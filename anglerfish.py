from __future__ import absolute_import
from __future__ import print_function
import argparse
import logging
from demux.demux import cluster_bc_matches, run_minimap2
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

    with open("adaptors.fasta", "w") as f:
        f.write(ss.get_fastastring())

    retcode = run_minimap2(args.in_fastq, "adaptors.fasta", "out.paf")

    for sample, adaptor in ss.samplesheet:
        cluster_bc_matches(args.in_fastq, args.out_fastq+"_"+sample+".fastq.gz", "out.paf", adaptor, args.max_distance, args.debug, args.count)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tools to demux I7 and I5 barcodes when sequenced by single-molecules")
    parser.add_argument('--in_fastq', '-i', required=True, help="Input ONT fastq file")
    parser.add_argument('--out_fastq', '-o', help="Output demux fastq file")
    parser.add_argument('--samplesheet', '-s', required=True, help="CSV formatted list of samples and barcodes")
    parser.add_argument('--count', '-c', action='store_true', help="Only do BC counting and not demuxing")
    parser.add_argument('--max-distance', '-m', default=1, type=int, help="Maximum edit distance for BC matching")
    parser.add_argument('--debug', '-d', action='store_true', help="Extra commandline output")
    args = parser.parse_args()
    # TODO: input validation
    #cluster_bc_matches(args.in_fastq, args.out_fastq, args.paf,args.i5_barcode,args.i7_barcode, args.max_distance, args.debug, args.count)
    run_demux(args)
