# Anglerfish

[![Anglerfish CI Status](https://github.com/remiolsen/anglerfish/workflows/Anglerfish/badge.svg)](https://github.com/remiolsen/anglerfish/actions)
[![PyPI](https://img.shields.io/pypi/v/bio-anglerfish)](https://pypi.python.org/pypi/bio-anglerfish/)
[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/anglerfish)](https://anaconda.org/bioconda/anglerfish)
[![Docker Container available](https://img.shields.io/docker/automated/remiolsen/anglerfish.svg)](https://hub.docker.com/r/remiolsen/anglerfish/)

## Introduction

Anglerfish is a tool designed to demultiplex Illumina libraries sequenced on Oxford Nanopore
flowcells. The primary purpose for this would be to do QC, i.e. to check pool balancing, assess
contamination, library insert sizes and so on.

For more information on how this can be used, please see this [poster](docs/AGBT_poster_20200214.pdf).

## Installation

### From PyPi

```
pip install bio-anglerfish
```

### From Bioconda

```
conda install -c bioconda anglerfish
```

### Install development version

```
pip install --upgrade --force-reinstall git+https://github.com/remiolsen/anglerfish.git
```

### For Arm64 processors (e.g. Apple M2)

Anglerfish depends on minimap2 which needs to be [compiled from source](https://github.com/lh3/minimap2?tab=readme-ov-file#installation) for Arm64 processors.
When minimap2 is compiled and available on $PATH, anglerfish can be installed with pip:

```shell
pip install bio-anglerfish
```

Additionaly, if Docker is your cup of tea, the Dockerfile supplied in this repository should also work on both arm64 and x86 processors.

## Source development

1. [Install miniconda](https://docs.conda.io/en/latest/miniconda.html).

2. Set up repo clone with editable install

For x86 processors (e.g. Intel/AMD):

```
git clone https://github.com/remiolsen/anglerfish.git
cd anglerfish
# Create a the anglerfish conda environment
conda env create -f environment.yml
# Install anglerfish
conda activate anglerfish-dev
pip install -e ".[dev]"
```

For Arm64 processors (e.g. Apple M1/M2). First [compile and install minimap2 manually](https://github.com/lh3/minimap2?tab=readme-ov-file#installation), then:

```
git clone https://github.com/remiolsen/anglerfish.git
cd anglerfish
# Create a the anglerfish conda environment (but remove minimap2)
conda env create -f <(grep -v minimap2 environment.yml)
# Install anglerfish
conda activate anglerfish-dev
pip install -e ".[dev]"
```

3. (Optional) Install pre-commit to prevent committing code that will fail linting

```
pre-commit install
```

4. (Optional) Enable automatic formatting in VS Code by creating `.vscode/settings.json` with:

```
{
  "editor.formatOnSave": true,
  "editor.defaultFormatter": "esbenp.prettier-vscode",
  "[python]": {
    "editor.defaultFormatter": "charliermarsh.ruff"
  },
  "prettier.configPath": "./pyproject.toml"
}
```

## Usage

Anglerfish requires two files to run.

- A basecalled FASTQ file from for instance Guppy (`/path/to/ONTreads.fastq.gz`)
- A samplesheet containing the sample names and indices expected to be found in the sequencing run. (`/path/to/samples.csv`)

Example of a samplesheet file:

```
P12864_201,truseq_dual,TAATGCGC-CAGGACGT,/path/to/ONTreads.fastq.gz
P12864_202,truseq_dual,TAATGCGC-GTACTGAC,/path/to/ONTreads.fastq.gz
P9712_101, truseq_dual,ATTACTCG-TATAGCCT,/path/to/ONTreads.fastq.gz
P9712_102, truseq_dual,ATTACTCG-ATAGAGGC,/path/to/ONTreads.fastq.gz
P9712_103, truseq_dual,ATTACTCG-CCTATCCT,/path/to/ONTreads.fastq.gz
P9712_104, truseq_dual,ATTACTCG-GGCTCTGA,/path/to/ONTreads.fastq.gz
P9712_105, truseq_dual,ATTACTCG-AGGCGAAG,/path/to/ONTreads.fastq.gz
P9712_106, truseq_dual,ATTACTCG-TAATCTTA,/path/to/ONTreads.fastq.gz
```

Or using single index (note samplesheet supports wildcard `*` use):

```
P12345_101,truseq,CAGGACGT,/path/to/*.fastq.gz
```

Then run:

```
anglerfish -s /path/to/samples.csv
```

### Options

#### Common

```
--out_fastq OUT_FASTQ, -o OUT_FASTQ
                      Analysis output folder (default: Current dir)
--samplesheet SAMPLESHEET, -s SAMPLESHEET
                      CSV formatted list of samples and barcodes
--threads THREADS, -t THREADS
                      Number of threads to use (default: 4)
--skip_demux, -c      Only do BC counting and not demuxing
--max-distance MAX_DISTANCE, -m MAX_DISTANCE
                       Manually set maximum edit distance for BC matching, automatically set this is set to either 1 or 2
--run_name RUN_NAME, -r RUN_NAME
                      Name of the run (default: anglerfish)
--debug, -d           Extra commandline output
--version, -v         Print version and quit

```

#### `--max-unknowns / -u`

Anglerfish will try to recover indices which are not specified in the samplesheet but follow the specified adaptor setup(s). This is analogous to `undetermined indices` as reported by Illumina demultiplexing. `--max-unknowns` will set the number of such indices reported.

#### `--lenient / -l`

This will consider both orientations of the I5 barcode and will use the reverse complement (of what was inputted in the samplesheet) only if significantly more reads were matched. This should be used with with extreme care, but the reason for this is that Anglerfish will try to guess which version of the Illumina samplesheet these indices were derived from. See this [guide](https://web.archive.org/web/20230602174828/https://knowledge.illumina.com/software/general/software-general-reference_material-list/000001800) for when i5 should be reverse complemented and not.

#### `--ont_barcodes / -n`

This is an ONT barcode aware mode. Which means each ONT barcode will be mapped and treated separately. A use case for this might be to put one Illumina pool per ONT barcode to spot potential index collisions you don't know of if you want to later make a pool of pools for sequencing in the same lane. This mode requires the fastq files to be placed in folders named `barcode01`, `barcode02`, etc. as is the default for MinKNOW (23.04). Example of such an anglerfish samplesheet:

```
P12345_101,truseq,CAGGACGT,/path/to/barcode01/*.fastq.gz
P54321_101,truseq,ATTACTCG,/path/to/barcode02/*.fastq.gz
```

### Output files

In folder `anglerfish_????_??_??_?????/`

- `*.fastq.gz` Demultiplexed reads (if any)
- `anglerfish_stats.txt` Barcode statistics from anglerfish run
- `anglerfish_stats.json` Machine readable anglerfish statistics

## Anglerfish Explore (Experimental)

`anglerfish-explore` is a command that aims to explore a sequencing pool without a given samplesheet and give hints on what adapter types are present, which index lenghts are used and whether there are any UMIs within the index sequence. The Anglerfish explore command is still under heavy development but can be triggered by running, e.g. for help text:

```shell
anglerfish-explore --help
```

## Credits

The Anglerfish code was written by [@remiolsen](https://github.com/remiolsen) but it would not exist without the contributions of [@FranBonath](https://github.com/FranBonath), [@taborsak](https://github.com/taborsak), [@ssjunnebo](https://github.com/ssjunnebo) and Carl Rubin.
Also, the [Anglerfish logo](docs/Anglerfish_logo.svg) was designed by [@FranBonath](https://github.com/FranBonath).

<p align="center">
  <img src="docs/Anglerfish_logo.svg">
</p>
