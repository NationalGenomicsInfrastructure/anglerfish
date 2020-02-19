# Anglerfish

[![Anglerfish CI Status](https://github.com/remiolsen/anglerfish/workflows/Anglerfish/badge.svg)](https://github.com/remiolsen/anglerfish/actions)
[![Docker Container available](https://img.shields.io/docker/automated/remiolsen/anglerfish.svg)](https://hub.docker.com/r/remiolsen/anglerfish/)

![Anglerfish](docs/Anglerfish_logo.svg)

## Introduction

Anglerfish is a tool designed to demultiplex Illumina libraries sequenced on Oxford Nanopore
flowcells. The primary purpose for this would be to do QC, i.e. to check pool balancing, assess
contamination, library insert sizes and so on.

## Installation

### Requirements

* Python3 (3.6)

Python modules:

* biopython v. 1.70
* python-levenshtein v. 0.12.0

Software:

* minimap2 v. 2.17
* fastqc v. 0.11.9 (optional)
* multiqc v. 1.8 (optional)

### Manually using pip

```
pip install --upgrade --force-reinstall git+https://github.com/remiolsen/anglerfish.git
```

### Manually with Conda

First [install miniconda](https://docs.conda.io/en/latest/miniconda.html), then:

```
git clone https://github.com/remiolsen/anglerfish.git
cd anglerfish
# Create a the anglerfish conda environment
conda env create -f environment.yml
# Install anglerfish
pip install -e .
```

### Docker

```
docker pull remiolsen/anglerfish:<version>
```

## Usage

Anglerfish requires two files to run.

  * A basecalled FASTQ file from for instance Guppy (`/path/to/ONTreads.fastq.gz`)
  * A samplesheet containing the sample names and indices expected to be found in the sequencing run. (`/path/to/samples.csv`)

Example of a samplesheet file:

```
P12864_201,truseq_dual,TAATGCGC,CAGGACGT
P12864_202,truseq_dual,TAATGCGC,GTACTGAC
P9712_101, truseq_dual,ATTACTCG,TATAGCCT
P9712_102, truseq_dual,ATTACTCG,ATAGAGGC
P9712_103, truseq_dual,ATTACTCG,CCTATCCT
P9712_104, truseq_dual,ATTACTCG,GGCTCTGA
P9712_105, truseq_dual,ATTACTCG,AGGCGAAG
P9712_106, truseq_dual,ATTACTCG,TAATCTTA
```

Or using single index:

```
P12345_101,truseq,CAGGACGT
```

Then run:

```
anglerfish.py -i /path/to/ONTreads.fastq.gz -o /path/to/samples.csv
```

### Optional

```
--out_fastq OUT_FASTQ, -o OUT_FASTQ
                      Analysis output folder (default: Current dir)
--samplesheet SAMPLESHEET, -s SAMPLESHEET
                      CSV formatted list of samples and barcodes
--threads THREADS, -t THREADS
                      Number of threads to use (default: 4)
--skip_demux, -c      Only do BC counting and not demuxing
--skip_fastqc, -f     After demuxing, skip running FastQC+MultiQC
--max-distance MAX_DISTANCE, -m MAX_DISTANCE
                      Maximum edit distance for BC matching (default: 2)
```

### Output files

In folder `anglerfish_????_??_??_?????/`

* `*.fastq.gz` Demultuplexed reads (if any)
* `anglerfish_stats.txt` Barcode statistics from anglerfish run
* `fastqc/` raw output from fastqc (if run)
* `multiqc/anglerfish_results_multiqc_report.html` Summary of demultiplexed reads
