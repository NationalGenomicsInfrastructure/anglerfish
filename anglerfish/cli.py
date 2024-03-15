import argparse
import os
from datetime import datetime as dt
from enum import Enum

import typer
from typing_extensions import Annotated

from .anglerfish import run_demux

app = typer.Typer(pretty_exceptions_show_locals=False)


class IndexOrientations(str, Enum):
    i7 = "i7"
    i5 = "i5"
    i7i5 = "i7+i5"
    default = "default"


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
            help="Manually set maximum edit distance for BC matching, automatically set this is set to either 1 or 2",
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
        int,
        typer.Option(
            "--lenient_factor",
            "-x",
            help="If lenient is set, this is the minimum factor of additional matches required to reverse complement the index",
        ),
    ] = 2,
    force_rc: Annotated[
        IndexOrientations,
        typer.Option(
            "--force_rc",
            "-p",
            help="Force reverse complementing the I5 and/or I7 indices. This will disregard lenient mode.",
        ),
    ] = IndexOrientations.default,
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
        bool, typer.Option("--version", "-v", help="Print version and quit")
    ] = False,
):
    """Run anglerfish demux. Now with emojis 💩✨"""
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
    utcnow = dt.utcnow()
    runname = utcnow.strftime(f"{args.run_name}_%Y_%m_%d_%H%M%S")
    assert os.path.exists(args.out_fastq)
    assert os.path.exists(args.samplesheet)
    args.out_fastq = os.path.join(os.path.abspath(args.out_fastq), runname)
    args.samplesheet = os.path.abspath(args.samplesheet)
    args.run_name = runname

    run_demux(args)
