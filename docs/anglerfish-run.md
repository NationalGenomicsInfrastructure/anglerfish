# Run demultiplexing, the standard method

This describes the the `anglerfish run` mode of running this tool, which expects an input samplesheet specifying the expected Illumina barcodes found in the sequenced ONT flowcell or barcode. Here we discuss some subjects in further details than the overview given in [Usage](#usage)

## Use cases

## Output

Example of file output from `angerfish run` with a single setup pool (as opposed the a [complex](#mixed-setup-pools) one) and without specifying an `--out_fastq` option thus generating a default name for the output folder.

```
anglerfish_run_YYYY_MM_DD_HHMMSS
├── anglerfish_dataframe.csv
├── anglerfish_stats.json
├── anglerfish_stats.txt
├── index_len(indexlength).fasta
└── index_len(indexlength).paf
```

The basic operation of Anglerfish is to map the input reads to a template of the adaptors `index_len(indexlength).fasta` using
minimap2 - output as an alignment file to `index_len(indexlength).paf`.

Anglerfish reports the stats of the run to a report called `anglerfish_stats.txt`, with the same number found in a machine readable JSON format (`anglerfish_stats.json`) Let's look at a few field from this report and number the lines:

```
01: Anglerfish v. 0.7.0 (run: anglerfish_2024_10_28_153312, 5c98ad62-784d-4b27-8dd4-a69bbfe553ac)
02: ===================
03: truseq_dual:
04: 105608	input_reads (100.00%)
05: 96593	reads aligning to adaptor sequences (91.46%)
06: 54785	aligned reads matching both I7 and I5 adaptor (56.72%)
07: 33658	aligned reads matching only I7 or I5 adaptor (34.85%)
08: 492	aligned reads matching multiple I7/I5 adaptor pairs (0.51%)
09: 7658	aligned reads with uncategorized alignments (7.93%)
```

- 03: Each adapter type will have their own section in the header
- 04: Any alignment from minimap given constraints of the [parameters](https://github.com/NationalGenomicsInfrastructure/anglerfish/blob/34ff1667d65281694e664bd48f53fa780f2075ce/anglerfish/demux/demux.py#L59) it's given
- 06: Reads matching the template (even partially) adapter1-insert-adapter2
- 08-09: Any [other matches](#uncategorized-alignments)

`anglerfish_dataframe.csv` is an attempt to summarize all index level stats (samplesheet samples and unknown indexes) into a
single "flat" table.

## Mixed setup pools

Anglerfish supports demultiplexing complex Illumina pools containing a mix of adapter setups, e.g. mixing samples with different index types like "truseq" types and "nextera" and samples with different lengths like 8+8bp and 6bp (single index).
See this [reference](https://web.archive.org/web/20231129095351/https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/AdapterSeq/Illumina_DNA/IlluminaUDIndexes.htm) from Illumina to get an idea of how these indexes might differ.

Example of such a samplesheet:

```
dual1,truseq_dual,TAATGCGC-CAGGACGT,/path/to/ONTreads.fastq.gz
dual2,truseq_dual,TAATGCGC-GTACTGAC,/path/to/ONTreads.fastq.gz
dual3,truseq_dual,ATTACTCG-TATAGCCT,/path/to/ONTreads.fastq.gz
single1,truseq,GAAACCCT,/path/to/ONTreads.fastq.gz
single2,truseq,CTGACTGA,/path/to/ONTreads.fastq.gz
single3,truseq,TCTCAGTG,/path/to/ONTreads.fastq.gz
```

The these are handled are, for each adapter-type and index length combination present, seperate minimap runs and read clustering is performed, then the results are aggregated in the report.
The path the fastq files supports glob'ing, e.g. you can specify multiple files like `/path/to/flowcell/fastq_passed/*.fastq.gz`

## Multiple ONT barcodes

## Unknown indexes

A list of indexes that do not match (within a set edit distance) of the indexes in the samplesheet will be listed in descending order at the bottom of the [report](#output).
These unknown matches are not clustered by sequence, such that each read error will get its' own entry therefore the list is truncated at `# samples in samplesheet` + 10. The column `closest_match` lists the samples(s) which have the shortest edit distance to this sequence.

The results might be distorted when the input fastq file(s) contain a [mixed adapter setup](#mixed-setup-pools). For the samplesheet given [above](#mixed-setup-pools) there unknown index list of both group of adaptor-index groups will be combined in the report. So the single-index group will contain hits to the "real", known dual index indices. E.g. "TAATGCGC" from "dual1" and "dual2" will be listed.

## Lenient mode

## Uncategorized alignments?
