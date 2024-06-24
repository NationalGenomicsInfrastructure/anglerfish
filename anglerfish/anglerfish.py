#!/usr/bin/env python
import glob
import gzip
import logging
import multiprocessing
import os
import sys
import uuid
from collections import Counter
from itertools import groupby

import Levenshtein as lev
import numpy as np
import pkg_resources

from .demux.adaptor import Adaptor
from .demux.demux import (
    Alignment,
    categorize_matches,
    cluster_matches,
    parse_reads_alns_from_paf,
    run_minimap2,
    write_demuxedfastq,
)
from .demux.report import AlignmentStat, Report, SampleStat
from .demux.samplesheet import SampleSheet

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("anglerfish")

MAX_PROCESSES = 64  # Ought to be enough for anybody

anglerfish_logo = r"""
     ___
   ( )  \ -..__
      _.|~”~~~”…_
    ^´           `>.
(+ (+ )             “<..<^(
  `´  ``´      ___       (
   \__..~      __(   _…_(
   \                /
    “--…_     _..~%´
         ```´´
"""


def run_demux(args):
    multiprocessing.set_start_method("spawn")

    if args.debug:
        log.setLevel(logging.DEBUG)
    run_uuid = str(uuid.uuid4())
    ss = SampleSheet(args.samplesheet, args.ont_barcodes)
    version = pkg_resources.get_distribution("bio-anglerfish").version
    report = Report(args.run_name, run_uuid, version)
    sys.stderr.write(anglerfish_logo)
    log.info(f" version {version}")
    log.info(f" arguments {vars(args)}")
    log.info(f" run uuid {run_uuid}")
    min_distance = ss.minimum_bc_distance()
    if args.max_distance is None:
        # Default: Set the maximum distance for barcode matching to 0, 1 or 2
        # depending on the smallest detected edit distance between indices in the samplesheet
        args.max_distance = min(min_distance - 1, 2)
        log.info(f"Using maximum edit distance of {args.max_distance}")
    if args.max_distance >= min_distance:
        log.error(
            f" The maximum allowed edit distance for barcode matching (={args.max_distance})"
            + f"is greater than the smallest detected edit distance between indices in samplesheet (={min_distance})"
            + ", which will result in ambiguous matches."
        )
        exit()
    log.debug(f"Samplesheet bc_dist == {min_distance}")
    if args.threads > MAX_PROCESSES:
        log.warning(
            f" Setting threads to {MAX_PROCESSES} as the maximum number of processes is {MAX_PROCESSES}"
        )
        args.threads = MAX_PROCESSES

    ## Sort the adaptors by type and size

    # Get a list of tuples with the adaptor name and ONT barcode
    adaptor_tuples: list[tuple[str, str]] = [
        (entry.adaptor.name, entry.ont_barcode) for entry in ss
    ]

    # Convert to set to enforce uniqueness
    adaptor_set: set[tuple[str, str]] = set(adaptor_tuples)

    # Create a dictionary with the adaptors as keys and an empty list as value
    adaptors_sorted: dict[tuple[str, str], list[tuple[str, Adaptor, str]]] = dict(
        [(i, []) for i in adaptor_set]
    )

    # Populate the dictionary values with sample-specific information
    """
        adaptors_sorted = {
             adaptor_name_str, ont_barcode_str ) : [
                (sample_name_str, Adaptor, fastq_str),
                (sample_name_str, Adaptor, fastq_str),
                ...
            ],
            ...
        } 
    """
    for entry in ss:
        adaptors_sorted[(entry.adaptor.name, entry.ont_barcode)].append(
            (entry.sample_name, entry.adaptor, os.path.abspath(entry.fastq))
        )
    if os.path.exists(args.out_fastq):
        raise FileExistsError(
            f"Output folder '{args.out_fastq}' already exists. Please remove it or specify another --run_name"
        )
    else:
        os.mkdir(args.out_fastq)
    out_fastqs = []
    for key, sample in adaptors_sorted.items():
        adaptor_name, ont_barcode = key
        fastq_path = sample[0][2]
        # If there are multiple ONT barcodes, we need to add the ONT barcode to the adaptor name
        if ont_barcode:
            adaptor_bc_name = f"{adaptor_name}_{ont_barcode}"
        else:
            adaptor_bc_name = adaptor_name
        fastq_files = glob.glob(fastq_path)

        # Align
        align_path = os.path.join(args.out_fastq, f"{adaptor_bc_name}.paf")
        adaptor_path = os.path.join(args.out_fastq, f"{adaptor_name}.fasta")
        with open(adaptor_path, "w") as f:
            f.write(ss.get_fastastring(adaptor_name))
        for fq in fastq_files:
            run_minimap2(fq, adaptor_path, align_path, args.threads)

        # Easy line count in input fastq files
        num_fq = 0
        for fq in fastq_files:
            with gzip.open(fq, "rb") as f:
                for i in f:
                    num_fq += 1
        num_fq = int(num_fq / 4)
        reads_alns: dict[str, list[Alignment]] = parse_reads_alns_from_paf(align_path)

        # Make stats
        log.info(f" Searching for adaptor hits in {adaptor_bc_name}")
        fragments, singletons, concats, unknowns = categorize_matches(
            adaptor_name + "_i5", adaptor_name + "_i7", reads_alns
        )
        stats = AlignmentStat(adaptor_bc_name)
        stats.compute_pafstats(num_fq, fragments, singletons, concats, unknowns)
        report.add_alignment_stat(stats)

        # Demux
        no_matches = []
        matches = []
        flipped_i7 = False
        flipped_i5 = False
        flips = {
            "original": {"i7_reversed": False, "i5_reversed": False},
            "i7": {"i7_reversed": True, "i5_reversed": False},
            "i5": {"i7_reversed": False, "i5_reversed": True},
            "i7+i5": {"i7_reversed": True, "i5_reversed": True},
        }
        if args.force_rc != "original":
            log.info(
                f" Force reverse complementing {args.force_rc} index for adaptor {adaptor_name}. Lenient mode is disabled"
            )
            no_matches, matches = cluster_matches(
                adaptors_sorted[key],
                fragments,
                args.max_distance,
                **flips[args.force_rc],
            )
            flipped_i7, flipped_i5 = flips[args.force_rc].values()
        elif args.lenient:  # Try reverse complementing the i5 and/or i7 indices and choose the best match
            flipped = {}
            results = []
            pool = multiprocessing.Pool(
                processes=4 if args.threads >= 4 else args.threads
            )
            results = []
            for flip, rev in flips.items():
                spawn = pool.apply_async(
                    cluster_matches,
                    args=(
                        adaptors_sorted[key],
                        fragments,
                        args.max_distance,
                        rev["i7_reversed"],
                        rev["i5_reversed"],
                    ),
                )
                results.append((spawn, flip))
            pool.close()
            pool.join()
            flipped = {result[1]: result[0].get() for result in results}

            best_flip = max(flipped, key=lambda k: len(flipped[k][1]))

            if (
                sorted([len(i[1]) for i in flipped.values()])[-1]
                == sorted([len(i[1]) for i in flipped.values()])[-2]
            ):
                log.warning(
                    " Lenient mode: Could not find any barcode reverse complements with unambiguously more matches. Using original index orientation for all adaptors. Please study the results carefully!"
                )
                no_matches, matches = flipped["original"]
            elif (
                best_flip != "None"
                and len(flipped[best_flip][1])
                > len(flipped["original"][1]) * args.lenient_factor
            ):
                log.info(
                    f" Lenient mode: Reverse complementing {best_flip} index for adaptor {adaptor_name} found at least {args.lenient_factor} times more matches"
                )
                no_matches, matches = flipped[best_flip]
                flipped_i7, flipped_i5 = flips[best_flip].values()
            else:
                log.info(
                    f" Lenient mode: using original index orientation for {adaptor_name}"
                )
                no_matches, matches = flipped["original"]
        else:
            no_matches, matches = cluster_matches(
                adaptors_sorted[key], fragments, args.max_distance
            )

        out_pool = []
        for k, v in groupby(sorted(matches, key=lambda x: x[3]), key=lambda y: y[3]):
            # To avoid collisions in fastq filenames, we add the ONT barcode to the sample name
            fq_prefix = k
            if ont_barcode:
                fq_prefix = ont_barcode + "-" + fq_prefix
            fq_name = os.path.join(args.out_fastq, fq_prefix + ".fastq.gz")
            out_fastqs.append(fq_name)
            sample_dict = {i[0]: [i] for i in v}

            # Find read lengths
            rlens = np.array([])
            for l, w in sample_dict.items():
                for i in w:
                    rlens = np.append(rlens, i[2] - i[1])
            rmean = np.round(np.mean(rlens), 2)
            rstd = np.round(np.std(rlens), 2)

            sample_stat = SampleStat(
                k,
                len(sample_dict.keys()),
                rmean,
                rstd,
                flipped_i7,
                flipped_i5,
                ont_barcode,
            )
            report.add_sample_stat(sample_stat)
            if not args.skip_demux:
                out_pool.append((sample_dict, fastq_path, fq_name))

        # Write demuxed fastq files
        pool = multiprocessing.Pool(processes=args.threads)
        results = []
        for out in out_pool:
            log.debug(f" Writing {out[2]}")
            spawn = pool.starmap_async(write_demuxedfastq, [out])
            results.append((spawn, out[2]))
        pool.close()
        pool.join()
        for result in results:
            log.debug(
                f" PID-{result[0].get()}: wrote {result[1]}, size {os.path.getsize(result[1])} bytes"
            )

        # Top unmatched indexes
        nomatch_count = Counter([x[3] for x in no_matches])
        if args.max_unknowns == 0:
            args.max_unknowns = len([sample for sample in ss]) + 10

        # We search for the closest sample in the samplesheet to the list of unknowns
        top_unknowns = []
        for i in nomatch_count.most_common(args.max_unknowns):
            sample_dists = [
                (
                    lev.distance(
                        i[0],
                        f"{x.adaptor.i7.index_seq}+{x.adaptor.i5.index_seq}".lower(),
                    ),
                    x.sample_name,
                )
                for x in ss
            ]
            closest_sample = min(sample_dists, key=lambda x: x[0])
            # If the distance is more than half the index length, we remove it
            if closest_sample[0] >= (len(i[0]) / 2) + 1:
                closest_sample = (closest_sample[0], None)
            else:
                # We might have two samples with the same distance
                all_min = [
                    x[1]
                    for x in sample_dists
                    if x[0] == closest_sample[0] and x[1] != closest_sample[1]
                ]
                # This list might be too long, so we truncate it
                if len(all_min) > 4:
                    all_min = all_min[:4]
                    all_min.append("...")
                if all_min:
                    closest_sample = (closest_sample[0], ";".join(all_min))

            top_unknowns.append([i[0], i[1], closest_sample[1]])
        report.add_unmatched_stat(top_unknowns, ont_barcode, adaptor_name)

    # Check if there were samples in the samplesheet without adaptor alignments and add them to report
    for entry in ss:
        if entry.sample_name not in [
            s.sample_name for s in [stat for stat in report.sample_stats]
        ]:
            sample_stat = SampleStat(entry.sample_name, 0, 0, 0, False, ont_barcode)
            report.add_sample_stat(sample_stat)

    report.write_report(args.out_fastq)
    report.write_json(args.out_fastq)
    report.write_dataframe(args.out_fastq, ss)
