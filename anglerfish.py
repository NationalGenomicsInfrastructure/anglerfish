from __future__ import absolute_import
from __future__ import print_function
import argparse
import logging
import sys
from itertools import groupby
from demux.demux import run_minimap2, parse_paf_lines, layout_matches, cluster_matches, write_demuxedfastq, run_fastqc
from demux.samplesheet import SampleSheet
from Bio.SeqIO.QualityIO import FastqGeneralIterator
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('anglerfish')


def run_demux(args):

    ss = SampleSheet(args.samplesheet)
    bc_dist = ss.minimum_bc_distance()
    if args.max_distance >= bc_dist:
        log.error(" Edit distance of barcodes in samplesheet are less than the minimum specified {}>={}".format(args.max_distance, bc_dist))
        exit()
    log.debug("Samplesheet bc_dist == {}".format(bc_dist))

    # TODO: split into several adaptor files and run minimap more than once
    with open("adaptors.fasta", "w") as f:
        f.write(ss.get_fastastring())
    retcode = run_minimap2(args.in_fastq, "adaptors.fasta", "out.paf", args.threads)

    # Sort the adaptors by type and size
    adaptors_t = [adaptor.name for sample, adaptor in ss]
    adaptor_set = set(adaptors_t)
    adaptors_sorted = dict([(i, []) for i in adaptor_set])
    for sample, adaptor in ss:
        adaptors_sorted[adaptor.name].append((sample, adaptor))

    paf_entries = parse_paf_lines("out.paf")
    out_fastqs = []
    for adaptor in adaptor_set:
        fragments, singletons, concats, unknowns = layout_matches(adaptor+"_i5",adaptor+"_i7",paf_entries)
        total = len(fragments)+len(singletons)+len(concats)+len(unknowns)
        log.info(" "+adaptor+" -- ")
        log.info(" {}\treads matching both I7 and I5 adaptor ({:.2f}%)".format(len(fragments), (len(fragments)/float(total)*100)))
        log.info(" {}\treads matching only I7 or I5 adaptor ({:.2f}%)".format(len(singletons), (len(singletons)/float(total)*100)))
        log.info(" {}\treads matching multiple I7/I5 adaptor pairs ({:.2f}%)".format(len(concats), (len(concats)/float(total)*100)))
        log.info(" {}\treads with uncategorized matches ({:.2f}%)".format(len(unknowns), (len(unknowns)/float(total)*100)))
        matches = cluster_matches(adaptors_sorted[adaptor], adaptor, fragments, args.max_distance)
        log.info("")
        log.info(" sample_name\t#reads")
        if not args.skip_demux:
            for k, v in groupby(sorted(matches,key=lambda x: x[3]), key=lambda y: y[3]):
                fq_name = k+".fastq.gz"
                out_fastqs.append(fq_name)
                sample_dict = {i[0]: [i] for i in v}
                log.info(" {}\t{}".format(k, len(sample_dict.keys())))
                write_demuxedfastq(sample_dict, args.in_fastq, fq_name)

    if not args.skip_fastqc and not args.skip_demux:
        fastqc = run_fastqc(out_fastqs, args.threads)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tools to demux I7 and I5 barcodes when sequenced by single-molecules')
    parser.add_argument('--in_fastq', '-i', required=True, help='Input ONT fastq file')
    parser.add_argument('--out_fastq', '-o', help='Output demux fastq file')
    parser.add_argument('--samplesheet', '-s', required=True, help='CSV formatted list of samples and barcodes')
    parser.add_argument('--threads', '-t', default=4, help='Number of threads to use')
    parser.add_argument('--skip_demux', '-c', action='store_true', help='Only do BC counting and not demuxing')
    parser.add_argument('--skip_fastqc', '-f', action='store_true', help='After demuxing, skip running FastQC+MultiQC')
    parser.add_argument('--max-distance', '-m', default=2, type=int, help='Maximum edit distance for BC matching')
    parser.add_argument('--debug', '-d', action='store_true', help='Extra commandline output')
    args = parser.parse_args()
    # TODO: input validation
    run_demux(args)
