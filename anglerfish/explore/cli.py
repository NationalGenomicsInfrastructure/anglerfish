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
    help="Number of threads specified by minimap2",
)
@click.option(
    "-n", "--no_overwrite", is_flag=True, help="Do not overwrite existing alignments"
)
@click.option(
    "-g",
    "--good_hit_threshold",
    default=0.9,
    type=float,
    help="Fraction of bases before and after index insert required to match perfectly for a hit to be considered a good hit. Default=0.9",
)
@click.option(
    "-i",
    "--insert_thres_low",
    default=4,
    type=int,
    help="Lower threshold for insert length, with value included",
)
@click.option(
    "-j",
    "--insert_thres_high",
    default=30,
    type=int,
    help="Upper threshold for insert length, with value included",
)
@click.option(
    "-B",
    "--minimap_b",
    default=4,
    type=int,
    help="Minimap2 -B parameter, mismatch penalty",
)
@click.option(
    "-m",
    "--min_hits_per_adaptor",
    default=50,
    type=int,
    help="Minimum number of good hits for an adaptor to be included in the analysis",
)
@click.option(
    "-u",
    "--umi_threshold",
    default=11,
    type=float,
    help="Minimum number of bases in insert to perform entropy calculation",
)
@click.option(
    "-k",
    "--kmer_length",
    default=2,
    type=int,
    help="Length of k-mers to use for entropy calculation",
)
def main(
    fastq,
    outdir,
    threads,
    no_overwrite,
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
        no_overwrite,
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
