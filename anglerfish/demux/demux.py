import glob
import io
import logging
import os
import re
import subprocess
from typing import cast

import Levenshtein as lev
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from anglerfish.demux.samplesheet import SampleSheetEntry

log = logging.getLogger("anglerfish")


def parse_cs(
    cs_string: str,
    index_seq: str,
    umi_before: int | None = 0,
    umi_after: int | None = 0,
) -> tuple[str, int]:
    """
    Given a cs string, an index sequence, and optional UMI lengths:

    - Parse the cs string to find the suspected index region of the read.
    - Return a tuple of the given index sequence and it's Levenshtein distance
        to the parsed index region of the read.
    """
    # Create pattern for a substitution from "n" in the adaptor to a base in the read
    n_subbed_pattern = re.compile(r"\*n([atcg])")
    # Concatenate all n-to-base substitutions to yield the sequence of the read spanning the mask
    bases_spanning_mask = "".join(re.findall(n_subbed_pattern, cs_string))
    # Trim away any UMIs
    if umi_before is not None and umi_before > 0:
        bases_spanning_index_mask = bases_spanning_mask[umi_before:]
    elif umi_after is not None and umi_after > 0:
        bases_spanning_index_mask = bases_spanning_mask[:-umi_after]
    else:
        bases_spanning_index_mask = bases_spanning_mask
    # Return the index and the Levenshtein distance between it and the presumed index region of the read
    return (
        bases_spanning_index_mask,
        lev.distance(index_seq.lower(), bases_spanning_index_mask),
    )


def run_minimap2(
    fastq_in: str,
    index_file: str,
    output_paf: str,
    threads: int,
    minimap_b: int = 1,
):
    """
    Runs Minimap2
    """
    cmd: list[str] = [
        "minimap2",
        "--cs",  # Output the cs tag (short)
        "-c",  # Output cigar string in .paf
        *["-A", "6"],  # Matching score
        *["-B", str(minimap_b)],  # Mismatch penalty
        *["-k", "10"],  # k-mer size
        *["-m", "8"],  # Minimal chaining score
        *["-t", str(threads)],  # Number of threads
        *["-w", "5"],  # Minimizer window size
        index_file,  # Target
        fastq_in,  # Query
    ]

    run_log = f"{output_paf}.log"
    with open(output_paf, "ab") as ofile, open(run_log, "ab") as log_file:
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log_file)
        subprocess.run("sort", stdin=p1.stdout, stdout=ofile, check=True)


class Alignment:
    """
    Class instantiated from one paf alignment entry.

    Retains the important values for later use.
    """

    def __init__(self, paf_line: str):
        paf_cols = paf_line.split()

        # Unpack cols to vars for type annotation
        self.read_name: str = paf_cols[0]
        self.read_len: int = int(paf_cols[1])  # read length
        self.read_start: int = int(paf_cols[2])  # start alignment on read
        self.read_end: int = int(paf_cols[3])  # end alignment on read
        self.read_strand: str = paf_cols[4]
        self.adapter_name: str = paf_cols[5]

        self.q: int = int(paf_cols[11])  # Q score

        self.cg: str = paf_cols[-2]  # cigar string
        self.cs: str = paf_cols[-1]  # cigar diff string


def map_reads_to_alns(
    paf_path: str, min_qual: int = 1, complex_identifier: bool = False
) -> dict[str, list[Alignment]]:
    """
    Parse a .paf file into a dict, mapping reads to their respective alignment objects.

    Outputs:

        reads_to_alns = {
            "read1": [
                aln_read1_adaptor1_i5,
                aln_read1_adaptor1_i7,
                aln_read1_adaptor2_i5,
                ...
            ],
            "read2": [
                aln_read1_adaptor1_i5,
                aln_read1_adaptor1_i7,
                aln_read1_adaptor2_i5,
                ...
            ],
            ...
        }

    complex_identifier = False (default)
        --> The keys will be on the form "{read_name}".

    complex_identifier = True
        --> The keys will be on the form "{read_name}_{i5_or_i7}_{strand_str}".
    """
    reads_to_alns: dict = {}
    with open(paf_path) as paf:
        for paf_line in paf:
            try:
                # Parse paf line into Alignment object
                aln = Alignment(paf_line)

                # Determine identifier
                i5_or_i7 = aln.adapter_name.split("_")[-1]
                strand_str = "positive" if aln.read_strand == "+" else "negative"
                key = (
                    f"{aln.read_name}_{i5_or_i7}_{strand_str}"
                    if complex_identifier
                    else aln.read_name
                )

            except IndexError:
                log.debug(f"Could not find all paf columns: {aln.read_name}")
                continue

            if aln.q < min_qual:
                log.debug(f"Low quality alignment: {aln.read_name}")
                continue

            if key in reads_to_alns.keys():
                reads_to_alns[key].append(aln)
            else:
                reads_to_alns[key] = [aln]

    return reads_to_alns


def categorize_matches(
    i5_name: str, i7_name: str, reads_to_alns: dict[str, list[Alignment]]
) -> tuple[
    dict[str, list[Alignment]],
    dict[str, list[Alignment]],
    dict[str, list[Alignment]],
    dict[str, list[Alignment]],
]:
    """
    Search the parsed paf alignments and layout possible Illumina library fragments,

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
    for read, alns in reads_to_alns.items():
        sorted_alns = []
        for i in range(len(alns) - 1):
            aln_i: Alignment = alns[i]
            aln_j: Alignment = alns[i + 1]
            if (
                aln_i.adapter_name != aln_j.adapter_name
                and (aln_i.adapter_name == i5_name and aln_j.adapter_name == i7_name)
                or (aln_j.adapter_name == i5_name and aln_i.adapter_name == i7_name)
            ):
                if aln_i in sorted_alns:
                    sorted_alns.append(aln_j)
                else:
                    sorted_alns.extend([aln_i, aln_j])
        if len(alns) == 1:
            singletons[read] = alns
        elif len(sorted_alns) == 2:
            fragments[read] = sorted(sorted_alns, key=lambda aln: aln.read_start)
        elif len(sorted_alns) > 2:
            concats[read] = sorted(sorted_alns, key=lambda aln: aln.read_start)
            log.debug(
                f"Concatenated fragment: {read}, found: {[(aln.adapter_name,aln.read_start) for aln in sorted_alns]}"
            )
        else:
            unknowns[read] = alns
            log.debug(
                f"Unknown fragment: {read}, found: {[(aln.adapter_name,aln.read_start) for aln in alns]}"
            )
        # TODO: add minimum insert size
    return (fragments, singletons, concats, unknowns)


def cluster_matches(
    entries: list[SampleSheetEntry],
    matches: dict[str, list[Alignment]],
    max_distance: int,
    i7_reversed: bool = False,
    i5_reversed: bool = False,
) -> tuple[list, list]:
    """Return the BED coordinates of unmatching and matching reads.

    Inputs:
    - matches: dict of reads to their respective alignments
    - entries: samplesheet entries for a given adaptor-barcode combination
    - max_distance: maximum allowed Levenshtein distance between index read and index sequence
    - i7_reversed: boolean indicating whether the i7 index is reversed
    - i5_reversed: boolean indicating whether the i5 index is reversed

    Outputs:
    (
        unmatched_bed: list of bed coordinates for unmatched reads
        matched_bed: list of bed coordinates for matched reads
    )
    """

    # Instantiate lists of bed coordinates to fill
    matched_bed = []
    unmatched_bed = []

    # Iterate over each read and it's alignments
    for read, alignments in matches.items():
        # Determine which alignment is i5 and which is i7
        if (
            alignments[0].adapter_name[-2:] == "i5"
            and alignments[1].adapter_name[-2:] == "i7"
        ):
            i5_aln = alignments[0]
            i7_aln = alignments[1]
        elif (
            alignments[1].adapter_name[-2:] == "i5"
            and alignments[0].adapter_name[-2:] == "i7"
        ):
            i5_aln = alignments[1]
            i7_aln = alignments[0]
        else:
            log.debug(f" {read} has no valid illumina fragment")
            continue

        entries_index_dists = []
        # Iterate over sample sheet entries and collect the distances between their index and the read's index
        for entry in entries:
            # Parse i5 index read
            if entry.adaptor.i5.index_seq is not None:
                if i5_reversed:
                    i5_seq = str(Seq(entry.adaptor.i5.index_seq).reverse_complement())
                else:
                    i5_seq = entry.adaptor.i5.index_seq

                i5_index_read, i5_index_dist = parse_cs(
                    i5_aln.cs,
                    i5_seq,
                    entry.adaptor.i5.len_umi_before_index,
                    entry.adaptor.i5.len_umi_after_index,
                )
            else:
                i5_index_read = ""
                i5_index_dist = 0

            # Parse i7 index read
            if entry.adaptor.i7.index_seq is not None:
                if i7_reversed:
                    i7_seq = str(Seq(entry.adaptor.i7.index_seq).reverse_complement())
                else:
                    i7_seq = entry.adaptor.i7.index_seq

                i7_index_read, i7_index_dist = parse_cs(
                    i7_aln.cs,
                    i7_seq,
                    entry.adaptor.i7.len_umi_before_index,
                    entry.adaptor.i7.len_umi_after_index,
                )
            else:
                i7_index_read = ""
                i7_index_dist = 0

            entries_index_dists.append(i5_index_dist + i7_index_dist)

        # Find the list idx of the minimum distance between the read's index and the indices of the samples in the sheet
        entries_min_index_dist_loc = min(
            range(len(entries_index_dists)), key=entries_index_dists.__getitem__
        )
        entry_min_index_dist = entries[entries_min_index_dist_loc]

        # If several samples in the sheet are equidistant from the read, skip the read
        if entries_index_dists.count(min(entries_index_dists)) > 1:
            continue

        # Get the coordinates of the read insert
        start_insert = min(i5_aln.read_end, i7_aln.read_end)
        end_insert = max(i7_aln.read_start, i5_aln.read_start)

        # If the read insert length is too short, skip the read
        if end_insert - start_insert < 10:
            continue

        # If the read is too far from the closest sample in the sheet
        if entries_index_dists[entries_min_index_dist_loc] > max_distance:
            # If the read has sensible lengths, add it to the unmatched beds
            if len(i7_index_read) + len(i5_index_read) == len(
                entry.adaptor.i7.index_seq or ""
            ) + len(entry.adaptor.i5.index_seq or ""):
                fi75 = "+".join(
                    [i for i in [i7_index_read, i5_index_read] if not i == ""]
                )
                unmatched_bed.append([read, start_insert, end_insert, fi75, "999", "."])
            # Otherwise, skip the read
            else:
                continue
        else:
            # Add the read to the matched beds
            matched_bed.append(
                [
                    read,
                    start_insert,
                    end_insert,
                    entry_min_index_dist.sample_name,
                    "999",
                    ".",
                ]
            )

    log.debug(
        f" Matched {len(matched_bed)} reads, unmatched {len(unmatched_bed)} reads"
    )

    return unmatched_bed, matched_bed


def write_demuxedfastq(
    beds: dict[str, list], fastq_in: os.PathLike, fastq_out: os.PathLike
) -> int:
    """
    Intended for multiprocessing
    Take a set of coordinates in bed format [[seq1, start, end, ..][seq2, ..]]
    from over a set of fastq entries in the input files and do extraction.

    Return: PID of the process
    """

    gz_buf = 131072
    fq_files = cast(list[str], glob.glob(fastq_in))
    assert len(fq_files) > 0, f"No fastq files found looking for {fastq_in}."
    for fq in fq_files:
        with subprocess.Popen(
            ["gzip", "-c", "-d", fq], stdout=subprocess.PIPE, bufsize=gz_buf
        ) as fzi:
            assert isinstance(fzi, subprocess.Popen)
            fi = io.TextIOWrapper(fzi.stdout, write_through=True)
            with open(fastq_out, "ab") as ofile:
                with subprocess.Popen(
                    ["gzip", "-c", "-f"],
                    stdin=subprocess.PIPE,
                    stdout=ofile,
                    bufsize=gz_buf,
                    close_fds=False,
                ) as oz:
                    assert isinstance(oz, subprocess.Popen)
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

    return os.getpid()
