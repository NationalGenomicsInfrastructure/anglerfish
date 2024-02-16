#!/usr/bin/env python
import argparse
import glob
import gzip
import logging
import os
import uuid
from collections import Counter
from datetime import datetime as dt
from itertools import groupby

import numpy as np
import pkg_resources

from .demux.demux import (
    cluster_matches,
    layout_matches,
    parse_paf_lines,
    run_minimap2,
    write_demuxedfastq,
)
from .demux.report import AlignmentStat, Report, SampleStat
from .demux.samplesheet import SampleSheet

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("anglerfish")


def run_demux(args):
    if args.debug:
        log.setLevel(logging.DEBUG)
    run_uuid = str(uuid.uuid4())
    os.mkdir(args.out_fastq)
    ss = SampleSheet(args.samplesheet, args.ont_barcodes)
    version = pkg_resources.get_distribution("bio-anglerfish").version
    report = Report(args.run_name, run_uuid, version)

    log.info(f" version {version}")
    log.info(f" arguments {vars(args)}")
    log.info(f" run uuid {run_uuid}")
    bc_dist = ss.minimum_bc_distance()
    if args.max_distance is None:
        if bc_dist > 1:
            args.max_distance = 2
        else:
            args.max_distance = 1
        log.info(f"Using maximum edit distance of {args.max_distance}")
    if args.max_distance >= bc_dist:
        log.error(
            f" Edit distance of barcodes in samplesheet are less than the minimum specified {args.max_distance}>={bc_dist}"
        )
        exit()
    log.debug(f"Samplesheet bc_dist == {bc_dist}")

    # Sort the adaptors by type and size
    adaptors_t = [(entry.adaptor.name, entry.ont_barcode) for entry in ss]
    adaptor_set = set(adaptors_t)
    adaptors_sorted = dict([(i, []) for i in adaptor_set])
    for entry in ss:
        adaptors_sorted[(entry.adaptor.name, entry.ont_barcode)].append(
            (entry.sample_name, entry.adaptor, os.path.abspath(entry.fastq))
        )

    out_fastqs = []
    for key, sample in adaptors_sorted.items():
        adaptor_name, ont_barcode = key
        fastq_path = sample[0][2]
        # If there are multiple ONT barcodes, we need to add the ONT barcode to the adaptor name
        adaptor_bc_name = adaptor_name
        if ont_barcode:
            adaptor_bc_name = adaptor_name + "_" + ont_barcode
        fastq_files = glob.glob(fastq_path)

        # Align
        aln_path = os.path.join(args.out_fastq, f"{adaptor_bc_name}.paf")
        adaptor_path = os.path.join(args.out_fastq, f"{adaptor_name}.fasta")
        with open(adaptor_path, "w") as f:
            f.write(ss.get_fastastring(adaptor_name))
        for fq in fastq_files:
            run_minimap2(fq, adaptor_path, aln_path, args.threads)

        # Easy line count in input fastq files
        num_fq = 0
        for fq in fastq_files:
            with gzip.open(fq, "rb") as f:
                for i in f:
                    num_fq += 1
        num_fq = int(num_fq / 4)
        paf_entries = parse_paf_lines(aln_path)

        # Make stats
        log.info(f" Searching for adaptor hits in {adaptor_bc_name}")
        fragments, singletons, concats, unknowns = layout_matches(
            adaptor_name + "_i5", adaptor_name + "_i7", paf_entries
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
            "i7": {"i7_reversed": True, "i5_reversed": False},
            "i5": {"i7_reversed": False, "i5_reversed": True},
            "i7+i5": {"i7_reversed": True, "i5_reversed": True},
        }
        if args.force_rc is not None:
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
        elif args.lenient:  # Try reverse complementing the I5 and/or i7 indices and choose the best match
            no_matches, matches = cluster_matches(
                adaptors_sorted[key], fragments, args.max_distance
            )
            flipped = {}
            for flip, rev in flips.items():
                rc_no_matches, rc_matches = cluster_matches(
                    adaptors_sorted[key], fragments, args.max_distance, **rev
                )
                flipped[flip] = (rc_matches, rc_no_matches, len(rc_matches))
            best_flip = max(flipped, key=lambda k: flipped[k][2])

            # There are no barcode flips with unambiguously more matches, so we abort
            if (
                sorted([i[2] for i in flipped.values()])[-1]
                == sorted([i[2] for i in flipped.values()])[-2]
            ):
                log.info(
                    "Could not find any barcode reverse complements with unambiguously more matches"
                )
            elif flipped[best_flip][2] > len(matches) * args.lenient_factor:
                log.info(
                    f" Reverse complementing {best_flip} index for adaptor {adaptor_name} found at least {args.lenient_factor} times more matches"
                )
                matches, no_matches, _ = flipped[best_flip]
                flipped_i7, flipped_i5 = flips[best_flip].values()
            else:
                log.info(f" Using original index orientation for {adaptor_name}")
        else:
            no_matches, matches = cluster_matches(
                adaptors_sorted[key], fragments, args.max_distance
            )

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
                write_demuxedfastq(sample_dict, fastq_path, fq_name)

        # Top unmatched indexes
        nomatch_count = Counter([x[3] for x in no_matches])
        if args.max_unknowns is None:
            args.max_unknowns = len([sample for sample in ss]) + 10
        report.add_unmatched_stat(
            nomatch_count.most_common(args.max_unknowns), ont_barcode, adaptor_name
        )

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

    if args.skip_fastqc:
        log.warning(
            " As of version 0.4.1, built in support for FastQC + MultiQC is removed. The '-f' flag is redundant."
        )


def anglerfish():
    parser = argparse.ArgumentParser(
        description="Tools to demux I7 and I5 barcodes when sequenced by single-molecules"
    )
    parser.add_argument(
        "--samplesheet",
        "-s",
        required=True,
        help="CSV formatted list of samples and barcodes",
    )
    parser.add_argument(
        "--out_fastq",
        "-o",
        default=".",
        help="Analysis output folder (default: Current dir)",
    )
    parser.add_argument(
        "--threads", "-t", default=4, help="Number of threads to use (default: 4)"
    )
    parser.add_argument(
        "--skip_demux",
        "-c",
        action="store_true",
        help="Only do BC counting and not demuxing",
    )
    parser.add_argument(
        "--skip_fastqc", "-f", action="store_true", help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--max-distance",
        "-m",
        type=int,
        help="Manually set maximum edit distance for BC matching, automatically set this is set to either 1 or 2",
    )
    parser.add_argument(
        "--max-unknowns",
        "-u",
        type=int,
        help="Maximum number of unknown indices to show in the output (default: length of samplesheet + 10)",
    )
    parser.add_argument(
        "--run_name",
        "-r",
        default="anglerfish",
        help="Name of the run (default: anglerfish)",
    )
    parser.add_argument(
        "--lenient",
        "-l",
        action="store_true",
        help="Will try reverse complementing the I5 and/or I7 indices and choose the best match.",
    )
    parser.add_argument(
        "--lenient_factor",
        "-x",
        default=4.0,
        type=float,
        help="If lenient is set, this is the minimum factor of additional matches required to reverse complement the index (default: 4.0)",
    )
    parser.add_argument(
        "--force_rc",
        "-p",
        choices=["i7", "i5", "i7+i5"],
        help="Force reverse complementing the I5 and/or I7 indices. This will disregard lenient mode.",
    )
    parser.add_argument(
        "--ont_barcodes",
        "-n",
        action="store_true",
        help="Will assume the samplesheet refers to a single ONT run prepped with a barcoding kit. And will treat each barcode separately",
    )
    parser.add_argument(
        "--debug", "-d", action="store_true", help="Extra commandline output"
    )
    parser.add_argument(
        "--version",
        "-v",
        action="version",
        help="Print version and quit",
        version=f'anglerfish {pkg_resources.get_distribution("bio-anglerfish").version}',
    )
    args = parser.parse_args()
    utcnow = dt.utcnow()
    runname = utcnow.strftime(f"{args.run_name}_%Y_%m_%d_%H%M%S")

    assert os.path.exists(args.out_fastq)
    assert os.path.exists(args.samplesheet)
    args.out_fastq = os.path.join(os.path.abspath(args.out_fastq), runname)
    args.samplesheet = os.path.abspath(args.samplesheet)
    args.run_name = runname
    run_demux(args)
