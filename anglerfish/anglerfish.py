#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
import argparse
import logging
import sys
import os
import json
import pkg_resources
import numpy as np
from datetime import datetime as dt
from itertools import groupby
from collections import Counter
from .demux.demux import run_minimap2, parse_paf_lines, layout_matches, cluster_matches, write_demuxedfastq
from .demux.samplesheet import SampleSheet
import gzip
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('anglerfish')


def run_demux(args):

    os.mkdir(args.out_fastq)
    ss = SampleSheet(args.samplesheet)
    version = pkg_resources.get_distribution("bio-anglerfish").version
    log.info(" version {}".format(version))
    log.info(" arguments {}".format(vars(args)))
    bc_dist = ss.minimum_bc_distance()
    if bc_dist == 0:
        log.error("There is one or more identical barcodes in the input samplesheet. Aborting!")
        exit()
    if args.max_distance == None:
        if bc_dist > 1:
            args.max_distance = 2
        else:
            args.max_distance = 1
        log.info("Using maximum edit distance of {}".format(args.max_distance))
    if args.max_distance >= bc_dist:
        log.error(" Edit distance of barcodes in samplesheet are less than the minimum specified {}>={}".format(args.max_distance, bc_dist))
        exit()
    log.debug("Samplesheet bc_dist == {}".format(bc_dist))

    # Sort the adaptors by type and size
    adaptors_t = [(adaptor.name, os.path.abspath(fastq)) for sample, adaptor, fastq in ss]
    adaptor_set = set(adaptors_t)
    adaptors_sorted = dict([(i, []) for i in adaptor_set])
    for sample, adaptor, fastq in ss:
        adaptors_sorted[(adaptor.name, os.path.abspath(fastq))].append((sample, adaptor))

    paf_stats = {}
    sample_stats = []
    unmatched_stats = []
    out_fastqs = []
    all_samples = []

    for key, sample in adaptors_sorted.items():

        adaptor_name, fastq_path = key
        sample0_name, adaptor0_object = sample[0]

        aln_name = "{}_{}".format(adaptor_name,os.path.basename(fastq_path).split('.')[0])
        assert aln_name not in paf_stats, "Input fastq files can not have the same filename"
        aln_path = os.path.join(args.out_fastq, "{}.paf".format(aln_name))
        adaptor_path = os.path.join(args.out_fastq,"{}.fasta".format(adaptor_name))
        with open(adaptor_path, "w") as f:
            f.write(ss.get_fastastring(adaptor_name))
        retcode = run_minimap2(fastq_path, adaptor_path, aln_path, args.threads)

        # Easy line count in input fastqfile
        fq_entries = 0
        with gzip.open(fastq_path, 'rb') as f:
            for i in f:
                fq_entries += 1
        fq_entries = int(fq_entries / 4)
        paf_entries = parse_paf_lines(aln_path)

        fragments, singletons, concats, unknowns = layout_matches(adaptor_name+"_i5",adaptor_name+"_i7",paf_entries)
        total = len(fragments)+len(singletons)+len(concats)+len(unknowns)

        paf_stats[aln_name] = {}
        paf_stats[aln_name]["input_reads"] = [fq_entries, 1.0]
        paf_stats[aln_name]["reads aligning to adaptor sequences"] = [total, total/float(fq_entries)]
        paf_stats[aln_name]["aligned reads matching both I7 and I5 adaptor"] = [len(fragments), len(fragments)/float(total)]
        paf_stats[aln_name]["aligned reads matching only I7 or I5 adaptor"] = [len(singletons), len(singletons)/float(total)]
        paf_stats[aln_name]["aligned reads matching multiple I7/I5 adaptor pairs"] = [len(concats), len(concats)/float(total)]
        paf_stats[aln_name]["aligned reads with uncategorized alignments"] = [len(unknowns), len(unknowns)/float(total)]

        no_matches, matches = cluster_matches(adaptors_sorted[key], adaptor_name, fragments, args.max_distance)

        aligned_samples = []
        for k, v in groupby(sorted(matches,key=lambda x: x[3]), key=lambda y: y[3]):
            aligned_samples.append(k)
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

            sample_stats.append("{}\t{}\t{}\t{}".format(k, len(sample_dict.keys()), rmean, rstd))
            if not args.skip_demux:
                write_demuxedfastq(sample_dict, fastq_path, fq_name)

        all_samples.extend(aligned_samples)
        # Check if there were samples in the samplesheet without adaptor alignments
        for ss_sample, ss_adaptor, ss_path in ss:
            if ss_adaptor.name == adaptor and ss_sample not in aligned_samples:
                sample_stats.append("{}\t0".format(ss_sample))

        # Top unmatched indexes
        nomatch_count = Counter([x[3] for x in no_matches])
        unmatched_stats.append(nomatch_count.most_common(10))

    header1 = ["sample_name","#reads","mean_read_len","std_read_len"]
    header2 = ["undetermined_index","count"]
    json_out = {"anglerfish_version":version,"paf_stats": [], "sample_stats": [], "undetermined": []}
    with open(os.path.join(args.out_fastq,"anglerfish_stats.txt"), "w") as f:
        f.write("Anglerfish v. "+version+"\n===================\n")
        for key, line in paf_stats.items():
            f.write(f"{key}:\n")
            for i,j in line.items():
                f.write(f"{j[0]}\t{i} ({j[1]*100:.2f}%)\n")
            json_out["paf_stats"].append(line)
        f.write("\n{}\n".format("\t".join(header1)))
        for sample in sample_stats:
            f.write(sample+"\n")
            json_out["sample_stats"].append(dict(zip(header1,sample.split("\t"))))
        f.write("\n{}\n".format("\t".join(header2)))
        for unmatch in unmatched_stats:
            for idx, mnum in unmatch:
                f.write("{}\t{}\n".format(idx, mnum))
                json_out["undetermined"].append(dict(zip(header2,[idx, mnum])))

    with open(os.path.join(args.out_fastq,"anglerfish_stats.json"), "w") as f:
        f.write(json.dumps(json_out,indent=2, sort_keys=True))
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
    parser.add_argument('--debug', '-d', action='store_true', help='Extra commandline output')
    args = parser.parse_args()
    utcnow = dt.utcnow()
    runname = utcnow.strftime("anglerfish_%Y_%m_%d_%H%M%S")

    assert os.path.exists(args.out_fastq)
    assert os.path.exists(args.samplesheet)
    args.out_fastq = os.path.join(os.path.abspath(args.out_fastq),runname)
    args.samplesheet = os.path.abspath(args.samplesheet)
    run_demux(args)
