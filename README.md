# Anglerfish

## Introduction

Anglerfish is a tool designed to demultiplex Illumina libraries sequenced on Oxford Nanopore
flowcells. The primary purpose for this would be to do QC, i.e. to check pool balancing, assess
contamination, library insert sizes and so on.

## Installation

Coming soon

## Usage

Anglerfish requires two files to run.

  * A basecalled FASTQ file from for instance Guppy
  * A samplesheet containing the sample names and indices expected to be found in the sequencing run.

Example of a samplesheet file:

```
Sample,    adaptors,   i5,      i7
P12864_201,truseq_dual,TAATGCGC,CAGGACGT
P12864_202,truseq_dual,TAATGCGC,GTACTGAC
P9712_101, truseq_dual,ATTACTCG,TATAGCCT
P9712_102, truseq_dual,ATTACTCG,ATAGAGGC
P9712_103, truseq_dual,ATTACTCG,CCTATCCT
P9712_104, truseq_dual,ATTACTCG,GGCTCTGA
P9712_105, truseq_dual,ATTACTCG,AGGCGAAG
P9712_106, truseq_dual,ATTACTCG,TAATCTTA
```
