# Introduction

These are the custom scripts/codes for the manuscript, "RNA is essential for PRC2 chromatin occupancy and function in human pluripotent stem cells" accepted in Nature Genetics 2020.

# ChIP-seq
### (1) Pipeline
See chipseq_pipeline.md.

### (2) Enrichment
Enrichment of genes was determined by DESeq2 using read counts in gene bodies. See chipseq_EZH2.md for example.

# rChIP-seq (RNase-dependent ChIP-seq)
### (1) Pipeline
same as ChIP-seq

### (2) Differential enrichment
Differential enrichment of genes between RNaseA treated and untreated was determined by DESeq2 using read counts in gene bodies. See rChIPseq_PRC2.md for example.

# RNA-seq
### (1) Differentially expressed genes
Differentially expressed genes were identified by DESeq2 using read counts in exons. See rnaseq_iPSC.md for example.
