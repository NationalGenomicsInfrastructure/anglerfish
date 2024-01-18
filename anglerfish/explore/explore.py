import logging
import os
import uuid

from anglerfish.demux.adaptor import load_adaptors
from anglerfish.demux.demux import run_minimap2

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("explore")


def run_explore(
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
    # Setup a run directory
    run_uuid = str(uuid.uuid4())
    try:
        os.mkdir(outdir)
    except FileExistsError:
        log.info(f"Output directory {outdir} already exists")
        if not no_overwrite:
            log.error(
                f"Output directory {outdir} already exists, please use --no_overwrite to continue"
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
        if os.path.exists(aln_path) and no_overwrite:
            log.info(f"Skipping {adaptor_name} as alignment already exists")
            continue
        adaptor_path = os.path.join(outdir, f"{adaptor_name}.fasta")
        with open(adaptor_path, "w") as f:
            f.write(adaptor.get_fastastring(insert_Ns=False))

        log.info(f"Aligning {adaptor_name}")
        run_minimap2(fastq, adaptor_path, aln_path, threads, minimap_b)
