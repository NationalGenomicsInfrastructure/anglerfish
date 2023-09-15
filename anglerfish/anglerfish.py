#!/usr/bin/env python
import argparse
import logging
import glob
import os
import pkg_resources
import numpy as np
import uuid
from datetime import datetime as dt
from itertools import groupby
from collections import Counter
from .demux.demux import run_minimap2, parse_paf_lines, layout_matches, cluster_matches, write_demuxedfastq
from .demux.samplesheet import SampleSheet
from .demux.report import Report, SampleStat, AlignmentStat
import gzip
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('anglerfish')


def run_demux(args):

    run_uuid = str(uuid.uuid4())
    os.mkdir(args.out_fastq)
    ss = SampleSheet(args.samplesheet)
    version = pkg_resources.get_distribution("bio-anglerfish").version
    report = Report(args.run_name, run_uuid, version)

    log.info(f" version {version}")
    log.info(f" arguments {vars(args)}")
    log.info(f" run uuid {run_uuid}")
    bc_dist = ss.minimum_bc_distance()
    if bc_dist == 0:
        log.error("There is one or more identical barcodes in the input samplesheet. Aborting!")
        exit()
    if args.max_distance == None:
        if bc_dist > 1:
            args.max_distance = 2
        else:
            args.max_distance = 1
        log.info(f"Using maximum edit distance of {args.max_distance}")
    if args.max_distance >= bc_dist:
        log.error(f" Edit distance of barcodes in samplesheet are less than the minimum specified {args.max_distance}>={bc_dist}")
        exit()
    log.debug(f"Samplesheet bc_dist == {bc_dist}")

    # Sort the adaptors by type and size
    adaptors_t = [(adaptor.name, os.path.abspath(fastq)) for _, adaptor, fastq in ss]
    adaptor_set = set(adaptors_t)
    adaptors_sorted = dict([(i, []) for i in adaptor_set])
    for sample, adaptor, fastq in ss:
        adaptors_sorted[(adaptor.name, os.path.abspath(fastq))].append((sample, adaptor))

    out_fastqs = []
    for key, sample in adaptors_sorted.items():

        adaptor_name, fastq_path = key
        fastq_files = glob.glob(fastq_path)

        aln_path = os.path.join(args.out_fastq, f"{adaptor_name}.paf")
        adaptor_path = os.path.join(args.out_fastq,f"{adaptor_name}.fasta")
        with open(adaptor_path, "w") as f:
            f.write(ss.get_fastastring(adaptor_name))
        for fq in fastq_files:
            retcode = run_minimap2(fq, adaptor_path, aln_path, args.threads)

        # Easy line count in input fastq files
        num_fq = 0
        for fq in fastq_files:
            with gzip.open(fq, 'rb') as f:
                for i in f:
                    num_fq  += 1
        num_fq  = int(num_fq  / 4)
        paf_entries = parse_paf_lines(aln_path)

        # Make stats
        fragments, singletons, concats, unknowns = layout_matches(adaptor_name+"_i5",adaptor_name+"_i7",paf_entries)
        stats = AlignmentStat(adaptor_name)
        stats.compute_pafstats(num_fq, fragments, singletons, concats, unknowns)
        report.add_alignment_stat(stats)

        no_matches, matches = cluster_matches(adaptors_sorted[key], fragments, args.max_distance)
        flipped = False
        if args.lenient:
            rc_no_matches, rc_matches = cluster_matches(adaptors_sorted[key], fragments, args.max_distance, i5_reversed=True)
            if len(rc_matches) * 0.2 > len(matches):
                log.info(f"Reversing I5 index for adaptor {adaptor_name} found at least 80% more matches")
                no_matches = rc_no_matches
                matches = rc_matches
                flipped = True

        for k, v in groupby(sorted(matches,key=lambda x: x[3]), key=lambda y: y[3]):
            fq_name = os.path.join(args.out_fastq, k+".fastq.gz")
            out_fastqs.append(fq_name)
            sample_dict = {i[0]: [i] for i in v}

            # Find read lengths
            rlens = np.array([])
            for l,w in sample_dict.items():
                for i in w:
                    rlens = np.append(rlens, i[2]-i[1])
            rmean = np.round(np.mean(rlens),2)
            rstd = np.round(np.std(rlens),2)

            sample_stat = SampleStat(k, len(sample_dict.keys()), rmean, rstd, flipped)
            report.add_sample_stat(sample_stat)
            if not args.skip_demux:
                write_demuxedfastq(sample_dict, fastq_path, fq_name)

        # Check if there were samples in the samplesheet without adaptor alignments and add them to report
        for ss_sample, _, _ in ss:
            if ss_sample not in [s.sample_name for s in [stat for stat in report.sample_stats]]:
                sample_stat = SampleStat(ss_sample, 0, 0, 0, False)
                report.add_sample_stat(sample_stat)

        # Top unmatched indexes
        nomatch_count = Counter([x[3] for x in no_matches])
        report.add_unmatched_stat(nomatch_count.most_common(args.max_unknowns))

    report.write_report(args.out_fastq)
    report.write_json(args.out_fastq)

    if args.skip_fastqc:
        log.warning(" As of version 0.4.1, built in support for FastQC + MultiQC is removed. The '-f' flag is redundant.")

def anglerfish():
    parser = argparse.ArgumentParser(description='Tools to demux I7 and I5 barcodes when sequenced by single-molecules')
    parser.add_argument('--samplesheet', '-s', required=True, help='CSV formatted list of samples and barcodes')
    parser.add_argument('--out_fastq', '-o', default='.', help='Analysis output folder (default: Current dir)')
    parser.add_argument('--threads', '-t', default=4, help='Number of threads to use (default: 4)')
    parser.add_argument('--skip_demux', '-c', action='store_true', help='Only do BC counting and not demuxing')
    parser.add_argument('--skip_fastqc', '-f', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--max-distance', '-m', type=int, help='Manually set maximum edit distance for BC matching, automatically set this is set to either 1 or 2')
    parser.add_argument('--max-unknowns', '-u', type=int, default=10, help='Maximum number of unknown indices to show in the output (default: 10)')
    parser.add_argument('--run_name', '-r', default='anglerfish', help='Name of the run (default: anglerfish)')
    parser.add_argument('--lenient', '-l', action='store_true', help='Will try reverse complementing the I5 index and choose the best match. USE WITH EXTREME CAUTION!')
    parser.add_argument('--debug', '-d', action='store_true', help='Extra commandline output')
    parser.add_argument('--version', '-v', action='version', help='Print version and quit', version=f'anglerfish {pkg_resources.get_distribution("bio-anglerfish").version}')
    args = parser.parse_args()
    utcnow = dt.utcnow()
    runname = utcnow.strftime(f"{args.run_name}_%Y_%m_%d_%H%M%S")

    assert os.path.exists(args.out_fastq)
    assert os.path.exists(args.samplesheet)
    args.out_fastq = os.path.join(os.path.abspath(args.out_fastq),runname)
    args.samplesheet = os.path.abspath(args.samplesheet)
    args.run_name = runname
    run_demux(args)
