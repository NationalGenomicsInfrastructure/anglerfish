from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import re
import logging
import Levenshtein as lev
import itertools
import subprocess
import io
from Bio.SeqIO.QualityIO import FastqGeneralIterator
logging.basicConfig(level=logging.INFO)


# TODO: Support for samplesheet to do real demuxing
# TODO: Build adapter fastas dynamically (based on idx length) and run minimap2
# TODO: un-hardcode this to enable different library types
trueseq_i5 = "I5"
trueseq_i7 = "I7"
trueseq_meta = {"I5": 70,
                "I7": 65
            }

def parse_cs(cs_string, index, max_distance):
    # Parses the CS string of a paf alignment and matches it to the given index using a max Levenshtein distance
    # TODO: Grab the alignment context and do Smith-Waterman,
    #       or do some clever stuff when parsing the cs string
    # PIPEDREAM: Do something big-brained with ONT squigglies
    log = logging.getLogger('demux')
    nt = re.compile("\*n([atcg])")
    nts = "".join(re.findall(nt, cs_string))
    #log.debug("Index search: {}".format(nts))

    # Allow for mismatches
    if lev.distance(index.lower(), nts) <= max_distance:
        return nts
    else:
        #log.debug(cs_string)
        return False



def cluster_bc_matches(in_fastq, out_fastq, paf, i5_barcode, i7_barcode, max_distance, debug, count):
    # Reads input ONT fastq file and clusters adaptor hits to find intact Illumina reads

    log = logging.getLogger('demux')
    if debug:
        log.info('Debug mode')
        log.setLevel(logging.DEBUG)

    # TODO: need to read more than 1 line
    with open(paf, "r") as p:
        oneln = p.readline()
        if not oneln.split()[23].startswith("cs:Z"):
            raise UserWarning("Input file was not a valid paf file or did not contain a cs:Z tag")
        if not oneln.split()[5] in trueseq_meta.keys():
            raise UserWarning("Input paf file does not aligned to a proper adapter sequences")

    raw_matches = {} # Raw matches from minimap2
    layout = {} # Sorted and index-annotated I7/I5 matches
    full_alns = 0 # DEBUG: all of adapter was mapped
    part_alns = 0 # DEBUG: partial mapping
    log.info("Parsing paf file")
    with open(paf, "r") as p:
        for line in p:
            aln = line.split()
            # TODO account for s2 tag. Or fiddle with minimap2 parameters
            try:
                entry = {"adapter": aln[5],
                         "rlen": int(aln[1]), # read length
                         "rstart": int(aln[2]), # start alignment on read
                         "rend": int(aln[3]), # end alignment on read
                         "strand": aln[4],
                         "cs": aln[23], # cs string
                         "iseq": None
                        }
            except IndexError:
                log.debug("Could not find all paf columns: {}".format(aln))

            if int(aln[3]) - int(aln[2]) != int(trueseq_meta[aln[5]]):
                part_alns += 1
            else:
                full_alns += 1

            if aln[0] in raw_matches.keys():
                raw_matches[aln[0]].append(entry)
            else:
                raw_matches[aln[0]] = [entry]

    log.debug("full alignments: {}, partial alignments {}".format(full_alns, part_alns))
    log.info("Searching for adaptor hits")
    index_match = 0
    no_index_match = 0
    for read, matches in raw_matches.items():
        ref_matches = []

        for match in matches:
            if match['adapter'] == trueseq_i5:
                found_i5 = parse_cs(match['cs'], i5_barcode, max_distance)
                if found_i5:
                    hit_i5 = match
                    hit_i5['iseq'] = found_i5

            elif match['adapter'] == trueseq_i7:
                found_i7 = parse_cs(match['cs'], i7_barcode, max_distance)
                if found_i7:
                    hit_i7 = match
                    hit_i7['iseq'] = found_i7

            ref_matches.append(match)

        ref_matches = sorted(ref_matches,key=lambda l:l['rstart'])
        layout[read] = ref_matches


    # Traverse layout and find intact Illumina reads.
    log.info("Finding Illumina reads")
    out_reads = {}; part_reads = {}
    for read, matches in layout.items():
        obed = []; pbed = []
        for match_i in range(1, len(matches), 2):
            m1 = matches[match_i-1]
            m2 = matches[match_i]
            bed_line = [read, m1['rend']+1, m2['rstart']-1, "insert_{}".format(match_i-1), "999", "."]

            if m1['iseq'] is not None and m2['iseq'] is not None:
                # We have consecutive I7+I5/I5+I7 matches
                # TODO: Check if these are actually not duplicates somehow, ie. I7+I7
                obed.append(bed_line)
                if m2['rstart']-m1['rend'] > 700:
                    log.debug("suspcious long read in: {}".format(read))
                    log.debug("{}".format(matches))
                if m2['rstart']-m1['rend'] < 30:
                    # TODO: Maybe drop these reads
                    log.debug("suspcious short read in: {}".format(read))
                    log.debug("{}".format(matches))

            elif m1['iseq'] is not None or m2['iseq'] is not None:
                pbed.append(bed_line)

        if len(obed) > 0:
            out_reads[read] = obed
        if len(pbed) > 0:
            part_reads[read] = pbed

    for read, bed in out_reads.items():
        log.debug("\t".join(str(x) for x in bed))
    log.info("demuxed reads {}".format(len(out_reads)))
    log.info("Reads with missing adaptor {}".format(len(part_reads)))

    if not count:
        write_demuxedfastq(out_reads, in_fastq, out_fastq)


def write_demuxedfastq(beds, fastq_in, fastq_out):
    # TODO: Can be optimized using pigz or rewritten using python threading
    gz_buf = 131072
    with subprocess.Popen(["gzip", "-c", "-d", fastq_in],
            stdout=subprocess.PIPE, bufsize=gz_buf) as fzi:
        fi = io.TextIOWrapper(fzi.stdout, write_through=True)
        with open(fastq_out, 'wb') as ofile:
            with subprocess.Popen(["gzip", "-c", "-f"],
                    stdin=subprocess.PIPE, stdout=ofile, bufsize=gz_buf, close_fds=False) as oz:

                for title, seq, qual in FastqGeneralIterator(fi):
                    new_title = title.split()
                    if new_title[0] not in beds.keys():
                        continue
                    outfqs = ""
                    for bed in beds[new_title[0]]:
                        new_title[0] += "_"+bed[3]
                        outfqs += "@{}\n".format(" ".join(new_title))
                        outfqs += "{}\n".format(seq[bed[1]:bed[2]])
                        outfqs += "+\n"
                        outfqs += "{}\n".format(qual[bed[1]:bed[2]])
                    oz.stdin.write(outfqs.encode('utf-8'))
