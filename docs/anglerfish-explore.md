# Explore demultiplexing, the experimental method

```{admonition} Note
  :kind: warning
  This mode of anglerfish work-in-progress and as such is not feature complete and might give sub-optimal results
```

This is a way of running anglerfish without using a samplesheet. The functionality is however quite different from `anglerfish run`
It will attempt, using read alignment, to find which adapter types are present in a given pool and their index length(s).

The basic usage of `anglerfish explore` expects two arguments, one for a single input fastq file `-f` and one for the output directory `-o`. Example shown for a pool where it finds `truseq_dual` adapters with 8+8 bp index length.

```
anglerfish explore -f fastq_passed/concat.fastq.gz -o explore
INFO:explore:Running anglerfish explore
INFO:explore:Run uuid 6b1c9b21-029f-4fc3-b28a-9214c776048a
INFO:explore:Aligning illumina_ud
INFO:explore:Aligning truseq
INFO:explore:Aligning truseq_dual
INFO:explore:Aligning truseq_umi
INFO:explore:Aligning nextera_legacy
INFO:explore:Aligning nextera_dual
INFO:explore:Parsing illumina_ud
INFO:explore:illumina_ud:i5 had 0 good hits.
INFO:explore:illumina_ud:i7 had 0 good hits.
INFO:explore:Adaptor illumina_ud is excluded from the analysis
INFO:explore:Parsing truseq
INFO:explore:truseq:i5 had 0 good hits.
INFO:explore:truseq:i7 had 5810 good hits.
INFO:explore:Adaptor truseq is excluded from the analysis
INFO:explore:Parsing truseq_dual
INFO:explore:truseq_dual:i5 had 779 good hits.
INFO:explore:truseq_dual:i7 had 5810 good hits.
INFO:explore:Adaptor truseq_dual is included in the analysis
INFO:explore:Parsing truseq_umi
INFO:explore:truseq_umi:i5 had 779 good hits.
INFO:explore:truseq_umi:i7 had 0 good hits.
INFO:explore:Adaptor truseq_umi is excluded from the analysis
INFO:explore:Parsing nextera_legacy
INFO:explore:nextera_legacy:i5 had 0 good hits.
INFO:explore:nextera_legacy:i7 had 0 good hits.
INFO:explore:Adaptor nextera_legacy is excluded from the analysis
INFO:explore:Parsing nextera_dual
INFO:explore:nextera_dual:i5 had 0 good hits.
INFO:explore:nextera_dual:i7 had 0 good hits.
INFO:explore:Adaptor nextera_dual is excluded from the analysis
INFO:explore:truseq_dual:i5 had 779 good hits with median insert length 8.0
INFO:explore:truseq_dual:i5 insert length histogram saved explore/truseq_dual_i5.hist.csv
INFO:explore:truseq_dual:i7 had 5810 good hits with median insert length 8.0
INFO:explore:truseq_dual:i7 insert length histogram saved explore/truseq_dual_i7.hist.csv
```
