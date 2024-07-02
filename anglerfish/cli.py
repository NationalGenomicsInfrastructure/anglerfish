import argparse
import datetime as dt
import os
from enum import Enum
from importlib.metadata import version as get_version
from typing import Optional

import typer
from typing_extensions import Annotated

from .anglerfish import run_demux
from .explore.explore import run_explore

app = typer.Typer(pretty_exceptions_show_locals=False)


class IndexOrientations(str, Enum):
    i7 = "i7"
    i5 = "i5"
    i7i5 = "i7+i5"
    original = "original"


def version_callback(value: bool):
    if value:
        print(f'anglerfish {get_version("bio-anglerfish")}')
        raise typer.Exit()


def deprecated_callback(value: bool):
    if value:
        raise typer.BadParameter(
            "Please use the 'anglerfish run -s' command to run anglerfish with a samplesheet. Running only 'anglerfish -s' is not supported as of version 0.7.0"
        )


@app.callback()
def main(
    version: Annotated[
        Optional[bool],
        typer.Option(
            "--version",
            "-v",
            help="Print version and quit",
            is_eager=True,
            callback=version_callback,
        ),
    ] = False,
    samplesheet: Annotated[
        Optional[str],
        typer.Option(
            "--samplesheet",
            "-s",
            hidden=True,
            is_eager=True,
            callback=deprecated_callback,
        ),
    ] = "",
):
    """
    Anglerfish is a tool designed to demultiplex Illumina libraries sequenced on Oxford Nanopore flowcells.
    The primary purpose for this would be to do QC, i.e. to check pool balancing, assess contamination, library insert sizes and so on.
    """
    if samplesheet:
        raise typer.BadParameter(
            "Please use the 'run' command to run anglerfish with a samplesheet. Running only 'anglerfish' is not supported as of version 0.7.0"
        )


@app.command()
def explore(
    fastq: Annotated[str, typer.Option("--fastq", "-f", help="Fastq file to align")],
    outdir: Annotated[str, typer.Option("--outdir", "-o", help="Output directory")],
    threads: Annotated[
        int,
        typer.Option(
            "--threads",
            "-t",
            help="Number of threads specified to minimap2",
        ),
    ] = 4,
    use_existing: Annotated[
        bool,
        typer.Option(
            "--use-existing",
            "-e",
            help="Use existing alignments if found in the specified output directory.",
        ),
    ] = False,
    good_hit_threshold: Annotated[
        float,
        typer.Option(
            "--good_hit_threshold",
            "-g",
            help="Fraction of adaptor bases immediately before and immediately after index insert required to match perfectly for a hit to be considered a good hit",
        ),
    ] = 0.9,
    insert_thres_low: Annotated[
        int,
        typer.Option(
            "--insert_thres_low",
            "-i",
            help="Lower threshold for index(+UMI) insert length, with value included.",
        ),
    ] = 4,
    insert_thres_high: Annotated[
        int,
        typer.Option(
            "--insert_thres_high",
            "-j",
            help="Upper threshold for index(+UMI) insert length, with value included.",
        ),
    ] = 30,
    minimap_b: Annotated[
        int,
        typer.Option(
            "--minimap_b",
            "-B",
            help="Minimap2 -B parameter, mismatch penalty.",
        ),
    ] = 4,
    min_hits_per_adaptor: Annotated[
        int,
        typer.Option(
            "--min_hits_per_adaptor",
            "-m",
            help="Minimum number of good hits for an adaptor to be included in the analysis.",
        ),
    ] = 50,
    umi_threshold: Annotated[
        float,
        typer.Option(
            "--umi_threshold",
            "-u",
            help="Minimum number of bases in insert to perform entropy calculation.",
        ),
    ] = 11,
    kmer_length: Annotated[
        int,
        typer.Option(
            "--kmer_length",
            "-k",
            help="Kmer length for entropy calculation.",
        ),
    ] = 2,
    version: Annotated[
        Optional[bool],
        typer.Option(
            "--version",
            "-v",
            help="Print version and quit",
            is_eager=True,
            callback=version_callback,
        ),
    ] = False,
):
    """This is an advanced samplesheet-free version of anglerfish."""
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


@app.command()
def run(
    samplesheet: Annotated[
        str,
        typer.Option(
            "--samplesheet", "-s", help="CSV formatted list of samples and barcodes"
        ),
    ],
    out_fastq: Annotated[
        str, typer.Option("--out_fastq", "-o", help="Analysis output folder")
    ] = ".",
    threads: Annotated[
        int, typer.Option("--threads", "-t", help="Number of threads to use")
    ] = 4,
    skip_demux: Annotated[
        bool,
        typer.Option("--skip_demux", "-c", help="Only do BC counting and not demuxing"),
    ] = False,
    max_distance: Annotated[
        int,
        typer.Option(
            "--max-distance",
            "-m",
            help="Manually set maximum allowed edit distance for index matching,"
            + "by default this is set to 0, 1 or 2 based on the minimum detected"
            + "index distance in the samplesheet.",
        ),
    ] = 2,
    max_unknowns: Annotated[
        int,
        typer.Option(
            "--max-unknowns",
            "-u",
            help="Maximum number of unknown indices to show in the output. default is length of samplesheet + 10",
        ),
    ] = 0,
    run_name: Annotated[
        str, typer.Option("--run_name", "-r", help="Run name")
    ] = "anglerfish",
    lenient: Annotated[
        bool,
        typer.Option(
            "--lenient",
            "-l",
            help="Will try reverse complementing the I5 and/or I7 indices and choose the best match.",
        ),
    ] = False,
    lenient_factor: Annotated[
        float,
        typer.Option(
            "--lenient_factor",
            "-x",
            help="If lenient is set, this is the minimum factor of additional matches required to reverse complement the index",
        ),
    ] = 4.0,
    force_rc: Annotated[
        IndexOrientations,
        typer.Option(
            "--force_rc",
            "-p",
            help="Force reverse complementing the I5 and/or I7 indices. If set to anything other than 'original' this will disregard lenient mode.",
        ),
    ] = IndexOrientations.original,
    ont_barcodes: Annotated[
        bool,
        typer.Option(
            "--ont_barcodes",
            "-n",
            help="Will assume the samplesheet refers to a single ONT run prepped with a barcoding kit. And will treat each barcode separately",
        ),
    ] = False,
    debug: Annotated[bool, typer.Option("--debug", "-d", help="Debug mode")] = False,
    version: Annotated[
        Optional[bool],
        typer.Option(
            "--version",
            "-v",
            help="Print version and quit",
            is_eager=True,
            callback=version_callback,
        ),
    ] = False,
):
    """Run anglerfish. This is the main command for anglerfish"""
    args = argparse.Namespace(
        samplesheet=samplesheet,
        out_fastq=out_fastq,
        threads=threads,
        skip_demux=skip_demux,
        max_distance=max_distance,
        max_unknowns=max_unknowns,
        run_name=run_name,
        lenient=lenient,
        lenient_factor=lenient_factor,
        force_rc=force_rc.value,
        ont_barcodes=ont_barcodes,
        debug=debug,
        version=version,
    )
    utcnow = dt.datetime.now(dt.timezone.utc)
    runname = utcnow.strftime(f"{args.run_name}_%Y_%m_%d_%H%M%S")
    if run_name != "anglerfish":
        runname = args.run_name

    assert os.path.exists(args.out_fastq), f"Output folder '{args.out_fastq}' not found"
    assert os.path.exists(
        args.samplesheet
    ), f"Samplesheet file '{args.samplesheet}' not found, please provide a valid path when using the --samplesheet option."
    args.out_fastq = os.path.join(os.path.abspath(args.out_fastq), runname)
    args.samplesheet = os.path.abspath(args.samplesheet)
    args.run_name = runname

    run_demux(args)
