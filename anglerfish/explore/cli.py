import click

from anglerfish.explore.explore import run_explore


@click.command()
@click.option("-f", "--fastq", required=True, help="Fastq file to align")
@click.option("-o", "--outdir", required=True, help="Output directory")
@click.option(
    "-t",
    "--threads",
    default=4,
    type=int,
    help="Number of threads specified to minimap2",
)
@click.option(
    "-e",
    "--use-existing",
    is_flag=True,
    help="Use existing alignments if found in the specified output directory.",
)
@click.option(
    "-g",
    "--good_hit_threshold",
    default=0.9,
    type=float,
    help="Fraction of adaptor bases immediately before and immediately after index insert required to match perfectly for a hit to be considered a good hit (default=0.9).",
)
@click.option(
    "-i",
    "--insert_thres_low",
    default=4,
    type=int,
    help="Lower threshold for index(+UMI) insert length, with value included (deafult=4).",
)
@click.option(
    "-j",
    "--insert_thres_high",
    default=30,
    type=int,
    help="Upper threshold for index(+UMI) insert length, with value included (default=30).",
)
@click.option(
    "-B",
    "--minimap_b",
    default=4,
    type=int,
    help="Minimap2 -B parameter, mismatch penalty (default=4).",
)
@click.option(
    "-m",
    "--min_hits_per_adaptor",
    default=50,
    type=int,
    help="Minimum number of good hits for an adaptor to be included in the analysis (default=50).",
)
@click.option(
    "-u",
    "--umi_threshold",
    default=11,
    type=float,
    help="Minimum number of bases in insert to perform entropy calculation (default=11).",
)
@click.option(
    "-k",
    "--kmer_length",
    default=2,
    type=int,
    help="Length of k-mers to use for entropy calculation (default=2).",
)
def main(
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
    run_explore(
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
    )


if __name__ == "__main__":
    main()
