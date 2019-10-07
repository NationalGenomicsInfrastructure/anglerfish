from __future__ import absolute_import
from __future__ import print_function
import argparse
import logging
import sys
from itertools import groupby
from demux.demux import run_minimap2, parse_paf_lines, layout_matches, cluster_matches, write_demuxedfastq
from demux.samplesheet import SampleSheet
from Bio.SeqIO.QualityIO import FastqGeneralIterator
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('anglerfish')


def run_demux(args):

    ss = SampleSheet(args.samplesheet)
    bc_dist = ss.minimum_bc_distance()
    if args.max_distance >= bc_dist:
        log.error("Edit distance of barcodes in samplesheet are less than the minimum specified {}>={}".format(args.max_distance, bc_dist))
        exit()
    log.debug("Samplesheet bc_dist == {}".format(bc_dist))

    # TODO: split into several adaptor files and run minimap more than once
    with open("adaptors.fasta", "w") as f:
        f.write(ss.get_fastastring())
    retcode = run_minimap2(args.in_fastq, "adaptors.fasta", "out.paf")

    # Sort the adaptors by type and size
    adaptors_t = [adaptor.name for sample, adaptor in ss]
    adaptor_set = set(adaptors_t)
    adaptors_sorted = dict([(i, []) for i in adaptor_set])
    for sample, adaptor in ss:
        adaptors_sorted[adaptor.name].append((sample, adaptor))

    paf_entries = parse_paf_lines("out.paf")
    for adaptor in adaptor_set:
        fragments, singletons, concats, unknowns = layout_matches(adaptor+"_i5",adaptor+"_i7",paf_entries)
        matches = cluster_matches(adaptors_sorted[adaptor], adaptor, fragments, args.max_distance)
        print(adaptor, len(fragments) ,len(singletons),len(concats), len(unknowns))
        for k, v in groupby(sorted(matches,key=lambda x: x[3]), key=lambda y: y[3]):
            sample_dict = {i[0]: [i] for i in v}
            write_demuxedfastq(sample_dict, args.in_fastq, k+".fastq.gz")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tools to demux I7 and I5 barcodes when sequenced by single-molecules")
    parser.add_argument('--in_fastq', '-i', required=True, help="Input ONT fastq file")
    parser.add_argument('--out_fastq', '-o', help="Output demux fastq file")
    parser.add_argument('--samplesheet', '-s', required=True, help="CSV formatted list of samples and barcodes")
    parser.add_argument('--count', '-c', action='store_true', help="Only do BC counting and not demuxing")
    parser.add_argument('--max-distance', '-m', default=2, type=int, help="Maximum edit distance for BC matching")
    parser.add_argument('--debug', '-d', action='store_true', help="Extra commandline output")
    args = parser.parse_args()
    # TODO: input validation
    run_demux(args)
