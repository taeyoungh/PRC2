RNA-seq analysis of iPSC
================
Taeyoung Hwang, PhD

### Library

``` r
library("readxl")
library(ggplot2)
library(cowplot)
library(reshape2)
library(DESeq2)
library(fdrtool)
```

# 1\. Import data

### Sample sheet

``` r
sampleSheet <- data.frame(read_excel("./RNAseq/RNAseq1_sampleSheet.xlsx"))
```

### featureCounts output

``` r
fCounts <- read.table("./RNAseq/RNAseq1_fCounts_table.txt", sep="\t", header=T, stringsAsFactors=F)
```

### Gene annotation

``` r
gencode <- read.table("./Annotation/gencodeV27_gene.txt", header=T, sep="\t", stringsAsFactors = F)
```

# 2\. Make data neat

### featureCounts

``` r
fCounts.annot <- fCounts[,1:2]
fCounts <- fCounts[,-c(1,2)]
rownames(fCounts) <- fCounts.annot$GeneID

fCounts <- fCounts[,match(sampleSheet$Sample,colnames(fCounts))]
```

# 3\. Gene expression

### TPM

``` r
tpm <- fCounts/(fCounts.annot$Length/1000)
tpm <- t(t(tpm)/(colSums(tpm)/10^6))
```

# 4\. Differential expression

### DESeq2

Define factors

``` r
tempSheet <- sampleSheet
tempSheet$CellLineSuppressed <- tempSheet$CellLine
tempSheet$CellLineSuppressed <- replace(tempSheet$CellLine, which(tempSheet$CellLine=="WT-A" | tempSheet$CellLine=="MT-A"), 1)
tempSheet$CellLineSuppressed <- replace(tempSheet$CellLineSuppressed, which(tempSheet$CellLine=="WT-B" | tempSheet$CellLine=="MT-B"), 2)

tempSheet$CellLineSuppressed <- factor(tempSheet$CellLineSuppressed, levels=c("1", "2"))
tempSheet$Condition <- factor(tempSheet$Condition, levels=c("WT", "MT"))
```

Run
DESeq2

``` r
dds<-DESeqDataSetFromMatrix(countData=fCounts, colData=tempSheet, design= ~ Condition + Condition:CellLineSuppressed)
dds<-DESeq(dds[apply(counts(dds),1, max)>=3,])
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
res <- results(dds, contrast = c("Condition","MT","WT"))
```

Distribution of
p-value

``` r
hist(res$pvalue, main="RNA-seq (iPSCs) with replicates", xlab="DESeq2 pvalue")
```

![](RNAseq_iPSC_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

### Table

``` r
deg <- data.frame(gencode[match(rownames(res), gencode$gene_id),], res)
deg <- deg[order(deg$pvalue),]
```

Define significant genes

``` r
deg$sig <- "Insig"
deg$sig[which(deg$padj<0.05)] <- "Sig"
table(deg$sig)
```

    ## 
    ## Insig   Sig 
    ## 28411  1918

``` r
deg$sigFC <- "Insig"
deg$sigFC[which(deg$padj<0.05)] <- "Sig"
deg$sigFC[which(deg$padj<0.05 & abs(deg$log2FoldChange)>=log2(2))] <- "SigFC"
table(deg$sigFC)
```

    ## 
    ## Insig   Sig SigFC 
    ## 28411  1451   467

Add TPM mean

``` r
idx <- match(deg$gene_id, rownames(fCounts))
deg$tpmMeanWT <- rowMeans(tpm[idx, subset(tempSheet, Condition=="WT")$Sample])
deg$tpmMeanMT <- rowMeans(tpm[idx, subset(tempSheet, Condition=="MT")$Sample])
deg$tpmMeanLog2FC <- log2(deg$tpmMeanMT / deg$tpmMeanWT)
```

# 5\. Ouput

``` r
write.table(deg, "./RNAseq/RNAseq_iPSC_deg.txt", sep="\t", row.names=F, col.names=T, quote=F)
```

``` r
rnaseq1.sampleSheet <- sampleSheet
rnaseq1.fCount <- fCounts
rnaseq1.fCount.annot <- fCounts.annot
rnaseq1.tpm <- tpm
rnaseq1.deg <- deg

save(list=c("rnaseq1.sampleSheet", "rnaseq1.fCount", "rnaseq1.fCount.annot", "rnaseq1.tpm", "rnaseq1.deg"), file="./RNAseq/RNAseq1_deg.Rdata")
```

``` r
sessionInfo()
```

    ## R version 3.3.3 (2017-03-06)
    ## Platform: x86_64-apple-darwin13.4.0 (64-bit)
    ## Running under: macOS  10.14.6
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] fdrtool_1.2.15             DESeq2_1.14.1             
    ##  [3] SummarizedExperiment_1.4.0 Biobase_2.34.0            
    ##  [5] GenomicRanges_1.26.4       GenomeInfoDb_1.10.3       
    ##  [7] IRanges_2.8.2              S4Vectors_0.12.2          
    ##  [9] BiocGenerics_0.20.0        reshape2_1.4.3            
    ## [11] cowplot_0.9.2              ggplot2_2.2.1             
    ## [13] readxl_1.0.0              
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bit64_0.9-7          splines_3.3.3        Formula_1.2-3       
    ##  [4] latticeExtra_0.6-28  blob_1.2.1           cellranger_1.1.0    
    ##  [7] yaml_2.2.1           pillar_1.4.3         RSQLite_2.2.0       
    ## [10] backports_1.1.2      lattice_0.20-35      digest_0.6.25       
    ## [13] RColorBrewer_1.1-2   XVector_0.14.1       checkmate_1.8.5     
    ## [16] colorspace_1.4-1     htmltools_0.4.0      Matrix_1.2-12       
    ## [19] plyr_1.8.4           XML_3.99-0.3         pkgconfig_2.0.3     
    ## [22] genefilter_1.56.0    zlibbioc_1.20.0      xtable_1.8-4        
    ## [25] scales_1.1.0         BiocParallel_1.8.2   htmlTable_1.13.3    
    ## [28] tibble_2.1.3         annotate_1.52.1      nnet_7.3-12         
    ## [31] lazyeval_0.2.2       survival_2.40-1      magrittr_1.5        
    ## [34] crayon_1.3.4         memoise_1.1.0        evaluate_0.14       
    ## [37] foreign_0.8-69       tools_3.3.3          data.table_1.10.4-3 
    ## [40] lifecycle_0.1.0      stringr_1.4.0        locfit_1.5-9.1      
    ## [43] munsell_0.5.0        cluster_2.0.6        AnnotationDbi_1.36.2
    ## [46] rlang_0.4.5          grid_3.3.3           RCurl_1.95-4.10     
    ## [49] rstudioapi_0.11      htmlwidgets_1.5.1    bitops_1.0-6        
    ## [52] base64enc_0.1-3      rmarkdown_2.1        gtable_0.3.0        
    ## [55] DBI_1.1.0            R6_2.4.1             gridExtra_2.3       
    ## [58] knitr_1.28           bit_1.1-15.2         Hmisc_4.1-1         
    ## [61] stringi_1.4.6        Rcpp_1.0.4.6         geneplotter_1.52.0  
    ## [64] vctrs_0.2.3          rpart_4.1-11         acepack_1.4.1       
    ## [67] xfun_0.12
