from __future__ import absolute_import
from __future__ import print_function
import argparse
import logging
from demux.demux import cluster_bc_matches
from Bio.SeqIO.QualityIO import FastqGeneralIterator
logging.basicConfig(level=logging.INFO)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tools to demux I7 and I5 barcodes when sequenced by single-molecules")
    parser.add_argument('--in_fastq', '-i', help="Input ONT fastq file")
    parser.add_argument('--out_fastq', '-o', help="Output demux fastq file")
    parser.add_argument('--paf', '-p', type=str, help="Minimap paf file")
    parser.add_argument('--i5-barcode', '-i5', type=str, help="i5 barcode sequence")
    parser.add_argument('--i7-barcode', '-i7', type=str, help="i7 barcode sequence")
    parser.add_argument('--count', '-c', action='store_true', help="Only do BC counting and not demuxing")
    parser.add_argument('--max-distance', '-m', default=1, type=int, help="Maximum edit distance for BC matching")
    parser.add_argument('--debug', '-d', action='store_true', help="Extra commandline output")
    args = parser.parse_args()
    cluster_bc_matches(args.in_fastq, args.out_fastq, args.paf,args.i5_barcode,args.i7_barcode, args.max_distance, args.debug, args.count)
