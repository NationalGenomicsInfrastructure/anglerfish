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
log = logging.getLogger('demux')



def parse_cs(cs_string, index, max_distance):
    """
    Parses the CS string of a paf alignment and matches it to the given index using a max Levenshtein distance
    TODO / idea: Do something big-brained with ONT squigglies
    """
    nt = re.compile("\*n([atcg])")
    nts = "".join(re.findall(nt, cs_string))

    # Allow for mismatches
    return nts, lev.distance(index.lower(), nts)


def run_minimap2(fastq_in, indexfile, output_paf, threads):
    """
    Runs Minimap2
    """
    cmd = [
        "minimap2",
        "--cs",
        "-m8",
        "-k", "10",
        "-w", "5",
        "-B1",
        "-A6",
        "--dual=no",
        "-c",
        "-t", str(threads),
        "-o", output_paf,
        indexfile,
        fastq_in
    ]

    proc = subprocess.run(cmd, check=True)
    return proc.returncode


def parse_paf_lines(paf, min_qual=10):
    """
    Read and parse one paf alignment lines.
    Returns a dict with the import values for later use
    """
    entries = {}
    with open(paf, "r") as paf:
        for paf_line in paf:
            aln = paf_line.split()
            try:
                # TODO: objectify this
                entry = {"adapter": aln[5],
                         "rlen": int(aln[1]), # read length
                         "rstart": int(aln[2]), # start alignment on read
                         "rend": int(aln[3]), # end alignment on read
                         "strand": aln[4],
                         "cs": aln[-1], # cs string
                         "q": int(aln[11]), # Q score
                         "iseq": None,
                         "sample": None
                        }
                read = aln[0]
            except IndexError:
                log.debug("Could not find all paf columns: {}".format(read))
                continue

            if entry["q"] < min_qual:
                log.debug("Low quality alignment: {}".format(read))
                continue

            if read in entries.keys():
                entries[read].append(entry)
            else:
                entries[read] = [entry]

    return entries


def layout_matches(i5_name, i7_name, paf_entries):
    """
    Search the parsed paf alignments and layout possible Illumina library fragments
    Returns dicts:
        - fragments. Reads with one I7 and one I5
        - singletons. Reads with that only match either I5 or I7 adaptors
        - concats. Concatenated fragments. Fragments with several alternating I7, I5 matches
        - unknowns. Any other reads
    """

    log.info(" Searching for adaptor hits")
    fragments = {}; singletons = {}; concats = {}; unknowns = {}
    for read, entry_list in paf_entries.items():
        sorted_entries = []
        for k in range(len(entry_list)-1):
            entry_i = entry_list[k]; entry_j = entry_list[k+1]
            if entry_i['adapter'] != entry_j['adapter'] and \
                (entry_i['adapter'] == i5_name and entry_j['adapter'] == i7_name) or \
                (entry_j['adapter'] == i5_name and entry_i['adapter'] == i7_name):
                if entry_i in sorted_entries:
                    sorted_entries.append(entry_j)
                else:
                    sorted_entries.extend([entry_i, entry_j])
        if len(entry_list) == 1:
            singletons[read] = entry_list
        elif len(sorted_entries) == 2:
            fragments[read] = sorted(sorted_entries,key=lambda l:l['rstart'])
        elif len(sorted_entries) > 2:
            concats[read] = sorted(sorted_entries,key=lambda l:l['rstart'])
        else:
            unknowns[read] = entry_list
        #TODO: add minimum insert size
    return (fragments, singletons, concats, unknowns)


def cluster_matches(sample_adaptor, adaptor_name, matches, max_distance):

    # Only illumina fragments
    matched = {}; matched_bed = []; unmatched_bed = []
    for read, alignments in matches.items():

        i5 = False
        i7 = False
        if alignments[0]['adapter'][-2:] == 'i5' and alignments[1]['adapter'][-2:] == 'i7':
            i5 = alignments[0]
            i7 = alignments[1]
        elif alignments[1]['adapter'][-2:] == 'i5' and alignments[0]['adapter'][-2:] == 'i7':
            i5 = alignments[1]
            i7 = alignments[0]
        else:
            log.debug(" Read has no valid illumina fragment")
            continue

        dists = []
        fi5 = ""; fi7 = ""
        for sample, adaptor in sample_adaptor:
            try:
                fi5, d1 = parse_cs(i5['cs'], adaptor.i5_index, max_distance)
            except AttributeError:
                d1 = 0 # presumably there it's single index
            fi7, d2 = parse_cs(i7['cs'], adaptor.i7_index, max_distance)
            dists.append(d1+d2)

        index_min = min(range(len(dists)), key=dists.__getitem__)
        # Test if two samples in the sheet is equidistant to the i5/i7
        if len([i for i, j in enumerate(dists) if j==dists[index_min]]) > 1:
            log.debug(" Ambiguous alignment, skipping")
            continue
        start_insert = min(i5['rend'],i7['rend'])
        end_insert = max(i7['rstart'],i5['rstart'])
        if end_insert - start_insert < 10:
            log.debug(" Erroneous / overlapping adaptor matches")
            continue
        if dists[index_min] > max_distance:
            log.debug(" No match")
            # Find only full length i7(+i5) adaptor combos. Basically a list of "known unknowns"
            if len(fi7) + len(fi5) == len(adaptor.i7_index or "") + len(adaptor.i5_index or ""):
                fi75 = "+".join([i for i in [fi7, fi5] if not i == ""])
                unmatched_bed.append([read, start_insert, end_insert, fi75, "999", "."])
            continue
        matched[read] = alignments
        matched_bed.append([read, start_insert, end_insert, sample_adaptor[index_min][0], "999", "."])
    return unmatched_bed, matched_bed



def write_demuxedfastq(beds, fastq_in, fastq_out):
    # Take a set of coordinates in bed format [[seq1, start, end, ..][seq2, ..]]
    # from over a set of fastq entries in the input files and do extraction.
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
