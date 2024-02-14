import logging
import os
import uuid

import pandas as pd

from anglerfish.demux.adaptor import load_adaptors
from anglerfish.demux.demux import parse_paf_lines, run_minimap2
from anglerfish.explore.entropy import calculate_relative_entropy

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("explore")


def run_explore(
    fastq,
    outdir,
    threads,
    use_existing,
    good_hit_threshold,
    insert_thres_low,
    insert_thres_high,
    minimap_b,
    min_hits_per_adaptor,
    umi_threshold,
    kmer_length,
):
    # Setup a run directory
    run_uuid = str(uuid.uuid4())
    try:
        os.mkdir(outdir)
    except FileExistsError:
        log.info(f"Output directory {outdir} already exists")
        if not use_existing:
            log.error(
                f"Output directory {outdir} already exists, please use --use_existing to continue"
            )
            exit(1)
        else:
            pass

    log.info("Running anglerfish explore")
    log.info(f"Run uuid {run_uuid}")

    adaptors = load_adaptors()
    alignments = []

    # Map all reads against all adaptors
    for adaptor in adaptors:
        adaptor_name = adaptor.name

        # Align
        aln_path = os.path.join(outdir, f"{adaptor_name}.paf")
        alignments.append((adaptor, aln_path))
        if os.path.exists(aln_path) and use_existing:
            log.info(f"Skipping {adaptor_name} as alignment already exists")
            continue
        adaptor_path = os.path.join(outdir, f"{adaptor_name}.fasta")
        with open(adaptor_path, "w") as f:
            f.write(adaptor.get_fastastring(insert_Ns=False))

        log.info(f"Aligning {adaptor_name}")
        run_minimap2(fastq, adaptor_path, aln_path, threads, minimap_b)

    # Parse alignments
    entries = {}
    adaptors_included = []
    for adaptor, aln_path in alignments:
        log.info(f"Parsing {adaptor.name}")
        aln_dict_with_lists = parse_paf_lines(aln_path, complex_identifier=True)

        # Choose only the highest scoring alignment for each combination of read, adaptor end and strand
        aln_dict = dict(
            [(k, max(v, key=lambda x: x["q"])) for k, v in aln_dict_with_lists.items()]
        )

        df = pd.DataFrame.from_dict(aln_dict, orient="index")
        nr_good_hits = {}

        # Match Insert Match = mim
        # The cs string filter is quite strict, requiring 10+ perfect match before insert and 10+ perfect match after insert
        # The cg string allows for mismatches within the matching strings but no insertions or deletions
        # All cg string matches are also cs string matches (subset) but not vice versa
        mim_re_cs = (
            r"^cs:Z::[1-9][0-9]*\+([a,c,t,g]*):[1-9][0-9]*$"  # Match Insert Match = mim
        )
        mim_re_cg = r"^cg:Z:([0-9]*)M([0-9]*)I([0-9]*)M$"
        df_mim = df[df.cs.str.match(mim_re_cs)]

        # Extract the match lengths
        match_col_df = df_mim.cg.str.extract(mim_re_cg).rename(
            {0: "match_1_len", 1: "insert_len", 2: "match_2_len"}, axis=1
        )
        match_col_df = match_col_df.astype(
            {
                "match_1_len": "int32",
                "insert_len": "int32",
                "match_2_len": "int32",
            }
        )

        df_mim.loc[match_col_df.index, match_col_df.columns] = match_col_df

        for adaptor_end_name, adaptor_end in zip(
            ["i5", "i7"], [adaptor.i5, adaptor.i7]
        ):
            if adaptor_end.has_index():
                # Alignment thresholds
                before_thres = round(
                    adaptor_end.len_before_index_region() * good_hit_threshold
                )
                after_thres = round(
                    adaptor_end.len_after_index_region() * good_hit_threshold
                )
                insert_thres_low = insert_thres_low
                insert_thres_high = insert_thres_high

                requirements = (
                    (df_mim["adapter"] == f"{adaptor.name}_{adaptor_end_name}")
                    & (df_mim["match_1_len"] >= (before_thres))
                    & (df_mim["insert_len"] >= insert_thres_low)
                    & (df_mim["insert_len"] <= insert_thres_high)
                    & (df_mim["match_2_len"] >= (after_thres))
                )
                df_good_hits = df_mim[requirements]

                median_insert_length = df_good_hits["insert_len"].median()
                insert_lengths = df_good_hits["insert_len"].value_counts()
            else:
                m_re_cs = r"^cs:Z::([1-9][0-9]*)$"
                df_good_hits = df[df.cg.str.match(m_re_cs)]
                match_col_df = df_good_hits.cg.str.extract(m_re_cs).rename(
                    {0: "match_1_len"}, axis=1
                )
                match_col_df = match_col_df.astype({"match_1_len": "int32"})

                df_good_hits.loc[
                    match_col_df.index, match_col_df.columns
                ] = match_col_df

                thres = round(
                    (
                        adaptor_end.len_before_index_region()
                        + adaptor_end.len_after_index_region()
                    )
                    * good_hit_threshold
                )
                df_good_hits = df_good_hits[df_good_hits["match_1_len"] >= thres]

                median_insert_length = None
                insert_lengths = None

            # Filter multiple hits per read and adaptor end
            df_good_hits = df_good_hits.sort_values(
                by=["read", "match_1_len"], ascending=False
            ).drop_duplicates(subset=["read", "adapter"], keep="first")

            if adaptor.name not in entries.keys():
                entries[adaptor.name] = {}
            entries[adaptor.name][adaptor_end_name] = df_good_hits

            nr_good_hits[adaptor_end_name] = len(df_good_hits)
            log.info(
                f"{adaptor.name}:{adaptor_end_name} had {len(df_good_hits)} good hits."
            )

        if min(nr_good_hits.values()) >= min_hits_per_adaptor:
            log.info(f"Adaptor {adaptor.name} is included in the analysis")
            adaptors_included.append(adaptor)
        else:
            log.info(f"Adaptor {adaptor.name} is excluded from the analysis")

    # Print insert length info for adaptor types included in the analysis
    for adaptor in adaptors_included:
        for adaptor_end_name, adaptor_end in zip(
            ["i5", "i7"], [adaptor.i5, adaptor.i7]
        ):
            df_good_hits = entries[adaptor.name][adaptor_end_name]
            if adaptor_end.has_index():
                median_insert_length = df_good_hits["insert_len"].median()
                if median_insert_length > umi_threshold:
                    # Calculate entropies here
                    entropies = calculate_relative_entropy(
                        df_good_hits, kmer_length, median_insert_length
                    )
                    entropy_file = os.path.join(
                        outdir, f"{adaptor.name}_{adaptor_end_name}.entropy.csv"
                    )
                    pd.Series(entropies).to_csv(entropy_file, float_format="%.2f")
                    log.info(
                        f"{adaptor.name}:{adaptor_end_name}, relative entropy for k={kmer_length}, over the index saved to {entropy_file}"
                    )
                insert_lengths = df_good_hits["insert_len"].value_counts()
                log.info(
                    f"{adaptor.name}:{adaptor_end_name} had {len(df_good_hits)} good hits with median insert length {median_insert_length}"
                )
                histogram_file = os.path.join(
                    outdir, f"{adaptor.name}_{adaptor_end_name}.hist.csv"
                )
                insert_lengths[sorted(insert_lengths.index)].to_csv(histogram_file)
                log.info(
                    f"{adaptor.name}:{adaptor_end_name} insert length histogram saved {histogram_file}"
                )
            else:
                median_insert_length = None
                insert_lengths = None
                log.info(
                    f"{adaptor.name}:{adaptor_end_name} had {len(df_good_hits)} good hits (no insert length since no index)"
                )
