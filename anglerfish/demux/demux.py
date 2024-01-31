import glob
import io
import logging
import os
import re
import subprocess

import Levenshtein as lev
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

log = logging.getLogger("anglerfish")


def parse_cs(cs_string, index, umi_before=0, umi_after=0):
    """
    Parses the CS string of a paf alignment and matches it to the given index using a max Levenshtein distance
    """
    nt = re.compile("\*n([atcg])")
    nts = "".join(re.findall(nt, cs_string))
    if umi_before > 0:
        nts = nts[umi_before:]
    if umi_after > 0:
        nts = nts[:-umi_after]
    # Allow for mismatches
    return nts, lev.distance(index.lower(), nts)


def run_minimap2(fastq_in, indexfile, output_paf, threads, minimap_b=1):
    """
    Runs Minimap2
    """
    cmd = [
        "minimap2",
        "--cs",
        "-m8",
        "-k",
        "10",
        "-w",
        "5",
        "-A",
        "6",
        "-B",
        str(minimap_b),
        "-c",
        "-t",
        str(threads),
        indexfile,
        fastq_in,
    ]

    run_log = f"{output_paf}.log"
    with open(output_paf, "ab") as ofile, open(run_log, "ab") as log_file:
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log_file)
        subprocess.run("sort", stdin=p1.stdout, stdout=ofile, check=True)


def parse_paf_lines(paf, min_qual=1, complex_identifier=False):
    """
    Read and parse one paf alignment lines.
    Returns a dict with the import values for later use
    If complex_identifier is True (default False), the keys will be on the form
    "{read}_{i5_or_i7}_{strand_str}"
    """
    entries = {}
    with open(paf) as paf:
        for paf_line in paf:
            aln = paf_line.split()
            try:
                # TODO: objectify this
                entry = {
                    "read": aln[0],
                    "adapter": aln[5],
                    "rlen": int(aln[1]),  # read length
                    "rstart": int(aln[2]),  # start alignment on read
                    "rend": int(aln[3]),  # end alignment on read
                    "strand": aln[4],
                    "cg": aln[-2],  # cigar string
                    "cs": aln[-1],  # cs string
                    "q": int(aln[11]),  # Q score
                    "iseq": None,
                    "sample": None,
                }
                read = entry["read"]
                if complex_identifier:
                    i5_or_i7 = entry["adapter"].split("_")[-1]
                    if entry["strand"] == "+":
                        strand_str = "positive"
                    else:
                        strand_str = "negative"
                    ix = f"{read}_{i5_or_i7}_{strand_str}"
                else:
                    ix = read
            except IndexError:
                log.debug(f"Could not find all paf columns: {read}")
                continue

            if entry["q"] < min_qual:
                log.debug(f"Low quality alignment: {read}")
                continue

            if ix in entries.keys():
                entries[ix].append(entry)
            else:
                entries[ix] = [entry]

    return entries


def layout_matches(i5_name, i7_name, paf_entries):
    """
    Search the parsed paf alignments and layout possible Illumina library fragments
    Returns dicts:
        - fragments. Reads with one I7 and one I5
        - singletons. Reads with that only match either I5 or I7 adaptors
        - concats. Concatenated fragments. Fragments with several alternating I7, I5 matches
        - unknowns. Any other reads, but usually i5-i5 or i7-i7 matches
    """

    fragments = {}
    singletons = {}
    concats = {}
    unknowns = {}
    for read, entry_list in paf_entries.items():
        sorted_entries = []
        for k in range(len(entry_list) - 1):
            entry_i = entry_list[k]
            entry_j = entry_list[k + 1]
            if (
                entry_i["adapter"] != entry_j["adapter"]
                and (entry_i["adapter"] == i5_name and entry_j["adapter"] == i7_name)
                or (entry_j["adapter"] == i5_name and entry_i["adapter"] == i7_name)
            ):
                if entry_i in sorted_entries:
                    sorted_entries.append(entry_j)
                else:
                    sorted_entries.extend([entry_i, entry_j])
        if len(entry_list) == 1:
            singletons[read] = entry_list
        elif len(sorted_entries) == 2:
            fragments[read] = sorted(sorted_entries, key=lambda l: l["rstart"])
        elif len(sorted_entries) > 2:
            concats[read] = sorted(sorted_entries, key=lambda l: l["rstart"])
            log.debug(
                f"Concatenated fragment: {read}, found: {[(i['adapter'],i['rstart']) for i in sorted_entries]}"
            )
        else:
            unknowns[read] = entry_list
            log.debug(
                f"Unknown fragment: {read}, found: {[(i['adapter'],i['rstart']) for i in entry_list]}"
            )
        # TODO: add minimum insert size
    return (fragments, singletons, concats, unknowns)


def cluster_matches(
    sample_adaptor, matches, max_distance, i7_reversed=False, i5_reversed=False
):
    # Only illumina fragments
    matched = {}
    matched_bed = []
    unmatched_bed = []
    for read, alignments in matches.items():
        i5 = False
        i7 = False
        if (
            alignments[0]["adapter"][-2:] == "i5"
            and alignments[1]["adapter"][-2:] == "i7"
        ):
            i5 = alignments[0]
            i7 = alignments[1]
        elif (
            alignments[1]["adapter"][-2:] == "i5"
            and alignments[0]["adapter"][-2:] == "i7"
        ):
            i5 = alignments[1]
            i7 = alignments[0]
        else:
            log.debug(" {read} has no valid illumina fragment")
            continue

        dists = []
        fi5 = ""
        fi7 = ""
        for _, adaptor, _ in sample_adaptor:
            try:
                i5_seq = adaptor.i5.index
                if i5_reversed and i5_seq is not None:
                    i5_seq = str(Seq(i5_seq).reverse_complement())
                fi5, d1 = parse_cs(
                    i5["cs"],
                    i5_seq,
                    adaptor.i5_umi_before,
                    adaptor.i5_umi_after,
                )
            except AttributeError:
                d1 = 0  # presumably it's single index, so no i5

            i7_seq = adaptor.i7.index
            if i7_reversed and i7_seq is not None:
                i7_seq = str(Seq(i7_seq).reverse_complement())
            fi7, d2 = parse_cs(
                i7["cs"],
                i7_seq,
                adaptor.i7_umi_before,
                adaptor.i7_umi_after,
            )
            dists.append(d1 + d2)

        index_min = min(range(len(dists)), key=dists.__getitem__)
        # Test if two samples in the sheet is equidistant to the i5/i7
        if len([i for i, j in enumerate(dists) if j == dists[index_min]]) > 1:
            continue
        start_insert = min(i5["rend"], i7["rend"])
        end_insert = max(i7["rstart"], i5["rstart"])
        if end_insert - start_insert < 10:
            continue
        if dists[index_min] > max_distance:
            # Find only full length i7(+i5) adaptor combos. Basically a list of "known unknowns"
            if len(fi7) + len(fi5) == len(adaptor.i7.index or "") + len(
                adaptor.i5.index or ""
            ):
                fi75 = "+".join([i for i in [fi7, fi5] if not i == ""])
                unmatched_bed.append([read, start_insert, end_insert, fi75, "999", "."])
            continue
        matched[read] = alignments
        matched_bed.append(
            [read, start_insert, end_insert, sample_adaptor[index_min][0], "999", "."]
        )
    log.debug(f" Matched {len(matched)} reads, unmatched {len(unmatched_bed)} reads")
    return unmatched_bed, matched_bed


def write_demuxedfastq(beds, fastq_in, fastq_out):
    """
    Take a set of coordinates in bed format [[seq1, start, end, ..][seq2, ..]]
    from over a set of fastq entries in the input files and do extraction.
    """
    gz_buf = 131072
    fq_files = glob.glob(fastq_in)
    for fq in fq_files:
        with subprocess.Popen(
            ["gzip", "-c", "-d", fq], stdout=subprocess.PIPE, bufsize=gz_buf
        ) as fzi:
            fi = io.TextIOWrapper(fzi.stdout, write_through=True)
            with open(fastq_out, "ab") as ofile:
                with subprocess.Popen(
                    ["gzip", "-c", "-f"],
                    stdin=subprocess.PIPE,
                    stdout=ofile,
                    bufsize=gz_buf,
                    close_fds=False,
                ) as oz:
                    for title, seq, qual in FastqGeneralIterator(fi):
                        new_title = title.split()
                        if new_title[0] not in beds.keys():
                            continue
                        outfqs = ""
                        for bed in beds[new_title[0]]:
                            new_title[0] += "_" + bed[3]
                            outfqs += "@{}\n".format(" ".join(new_title))
                            outfqs += f"{seq[bed[1] : bed[2]]}\n"
                            outfqs += "+\n"
                            outfqs += f"{qual[bed[1] : bed[2]]}\n"
                        oz.stdin.write(outfqs.encode("utf-8"))
        log.debug(f" Wrote {fastq_out}, size: {os.path.getsize(fastq_out)} bytes")
