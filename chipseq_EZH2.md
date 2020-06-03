ChIP-seq analysis - EZH2 enrichment
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
sampleSheet <- data.frame(read_excel("./ChIPseq/ChIPseq_summary.xlsx"))
```

    ## Warning in strptime(x, format, tz = tz): unknown timezone 'zone/tz/2019c.1.0/
    ## zoneinfo/America/Denver'

``` r
sampleSheet <- subset(sampleSheet, Experiment =="ChIPseq" & Library %in% c("Lysate_Input", "IP_EZH2") & Treatment=="None")
```

### featureCounts output

``` r
ip.fCounts.table <- read.table("./ChIPseq/ChIPseq_EZH2_fCounts.txt", header=T, sep="\t", stringsAsFactors = F, check.names = F)
input.fCounts.table <- read.table("./ChIPseq/ChIPseq_Input_fCounts.txt", header=T, sep="\t", stringsAsFactors = F, check.names = F)
```

### Gene annotation

``` r
gencode <- read.table("./Annotation/gencodeV27_gene.txt", header=T, sep="\t", stringsAsFactors = F)
```

# 2\. Make data neat

### featureCounts

IP

``` r
ip.fCounts.annot <- ip.fCounts.table[, c(1,2)]
ip.fCounts.count <- ip.fCounts.table[, -c(1,2)]
rownames(ip.fCounts.count) <- ip.fCounts.annot$GeneID
```

Input

``` r
input.fCounts.annot <- input.fCounts.table[, c(1,2)]
input.fCounts.count <- input.fCounts.table[, -c(1,2)]
rownames(input.fCounts.count) <- input.fCounts.annot$GeneID
```

Merge IP and Input

``` r
which(ip.fCounts.annot$GeneID != input.fCounts.annot$GeneID)
```

    ## integer(0)

``` r
which(ip.fCounts.annot$Length != input.fCounts.annot$Length)
```

    ## integer(0)

``` r
fCounts.annot <- ip.fCounts.annot
colnames(fCounts.annot) <- c("gene_id", "size")
fCounts.count <- cbind(ip.fCounts.count, input.fCounts.count)
```

### Pooled version of featureCounts

Prepare pooled
name

``` r
sampleSheet$CellLine <- factor(sampleSheet$CellLine, levels=c("WT_A", "WT_B", "MT_A", "MT_B"), labels=c("WTa", "WTb", "MTa", "MTb"))
sampleSheet$Pulldown <- sapply(strsplit(sampleSheet$Library, split="_"), "[[", 2)
sampleSheet$PooledName <- paste(sampleSheet$CellLine, sampleSheet$Pulldown, sep="_")
```

Perform
pooling

``` r
pooled.fCounts.count <- sapply(split(sampleSheet, sampleSheet$PooledName), function(x) return(rowSums(fCounts.count[, x$Sample, drop=F])))
```

### Pooled sheet

``` r
pooledSheet <- data.frame(CellLine=rep(c("WTa", "WTb", "MTa", "MTb"), each=2), Condition=rep(c("WT", "MT"), each=4), Pulldown=rep(c("EZH2", "Input"), 4), stringsAsFactors = F)
pooledSheet$PooledName <- paste(pooledSheet$CellLine, pooledSheet$Pulldown, sep="_")
```

Add dedupped read number for
normalization

``` r
temp <- dcast(sampleSheet, PooledName~., value.var = "DeduppedNum", fun.aggregate = sum)
pooledSheet$DeduppedReadNum <- temp[match(pooledSheet$PooledName, temp$PooledName), 2]
```

``` r
colnames(pooledSheet) <- c("CellLine", "Condition", "Library", "Sample", "DeduppedReadNum")
```

# 3\. Enrichment

### functions

Function for enrichment calculation using DESeq2 for a given
design

``` r
runDESeq2ForChIPseq <- function(sampleSheet, countMatrix, designFormula, normVector, count_th_fun="max", count_th_arg=3) {

  tempCount <- countMatrix[, match(sampleSheet$Sample, colnames(countMatrix))]
  rownames(sampleSheet) <- sampleSheet$Sample
  
  deseq2.obj <- DESeqDataSetFromMatrix(countData=tempCount, colData=sampleSheet , design=formula(designFormula))
  sizeFactors(deseq2.obj) <- normVector/min(normVector)
  deseq2.obj <- DESeq(deseq2.obj[apply(counts(deseq2.obj),1,get(count_th_fun))>=count_th_arg,])
  output <- results(deseq2.obj, contrast=c("Library", "Pulldown","Input")) 
  return(output)
}
```

design should be a form of " ~ Covariates + Library“.  
sampleSheet should have a column of”Library" consisting of “Input” and
“Pulldown”.  
sampleSheet should have a column of “Sample” consisting of unique
names.  
\- Normalization is done with a given normalization vector.  
\- Enrichment is calculated by contrasting pulldown vs input.

Function for p-value calculation based on empirical estimation of
variation

``` r
assignEmpPvalue <- function(input, main="Input") {
  par(mfrow=c(1,4))
  
  hist(input$pvalue, main=main, xlab="DESeq2 pvalue")
  hist(input$stat, main="Histogram", xlab="DESeq2 statistic")
  
  temp <- fdrtool(input$stat, statistic= "normal", plot = F)
  qqnorm(input$stat, main="Q-Q normal", ylab="DESeq2 statistic Quantiles")
  qqline(input$stat, col="red", lwd=2); abline(h=0,v=0,lty=2)
  
  hist(temp$pval, main="DESeq2, followed by fdrtool", xlab="pvalue")
  
  input$pvalue=temp$pval
  input$padj=temp$qval
  return(input)
}
```

Function for converting DESeq2 output to table

``` r
deseq2table <- function(input, annotTable, pvalTh) {
  # Significant genes
  input$sig="Insig"
  input$sig[which(input$padj<pvalTh)]="Sig"
  
  input$sigFC="Insig"
  input$sigFC[which(input$padj<pvalTh)]="Sig"
  input$sigFC[which(input$padj<pvalTh & abs(input$log2FoldChange)>1)]="SigFC"
  print(table(input$sigFC))
  
  # data.frame  
  input <- as.data.frame(input)
  out <- cbind(annotTable[match(rownames(input), annotTable$gene_id),], input)
  out <- out[order(out$pval), ]
  rownames(out) <- 1:nrow(out)
  return(out)
}
```

### WT enrichment

``` r
enrich <- list()
```

Define
factors

``` r
tempSheet <- subset(pooledSheet, Condition=="WT" & Library %in% c("EZH2", "Input"))
tempSheet$CellLine <- factor(tempSheet$CellLine)
tempSheet$Library <- factor(tempSheet$Library, levels=c("Input", "EZH2"), labels=c("Input", "Pulldown"))
```

DEseq2

``` r
enrich[["WT"]] <- runDESeq2ForChIPseq(tempSheet, pooled.fCounts.count, " ~ CellLine+Library", normVector=tempSheet$DeduppedReadNum)
```

    ## converting counts to integer mode

    ## using pre-existing size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

Empirical p-value

``` r
enrich[["WT"]] <- assignEmpPvalue(enrich[["WT"]], main="EZH2: WT")
```

    ## Step 1... determine cutoff point
    ## Step 2... estimate parameters of null distribution and eta0
    ## Step 3... compute p-values and estimate empirical PDF/CDF
    ## Step 4... compute q-values and local fdr

![](ChIPseq_EZH2_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

Table

``` r
enrich[["WT"]] <- deseq2table(enrich[["WT"]], gencode, 0.05)
```

    ## 
    ## Insig   Sig SigFC 
    ## 32944  4314   293

### MT enrichment

Define
factors

``` r
tempSheet <- subset(pooledSheet, Condition=="MT" & Library %in% c("EZH2", "Input"))
tempSheet$CellLine <- factor(tempSheet$CellLine)
tempSheet$Library <- factor(tempSheet$Library, levels=c("Input", "EZH2"), labels=c("Input", "Pulldown"))
```

DEseq2

``` r
enrich[["MT"]] <- runDESeq2ForChIPseq(tempSheet, pooled.fCounts.count, " ~ CellLine+Library", normVector=tempSheet$DeduppedReadNum)
```

    ## converting counts to integer mode

    ## using pre-existing size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

Empirical p-value

``` r
enrich[["MT"]] <- assignEmpPvalue(enrich[["MT"]], main="EZH2: MT")
```

    ## Step 1... determine cutoff point
    ## Step 2... estimate parameters of null distribution and eta0
    ## Step 3... compute p-values and estimate empirical PDF/CDF
    ## Step 4... compute q-values and local fdr

![](ChIPseq_EZH2_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

Table

``` r
enrich[["MT"]] <- deseq2table(enrich[["MT"]], gencode, 0.05)
```

    ## 
    ## Insig   Sig SigFC 
    ## 33121  4505    68

### Differential enrichment

Define factors

``` r
tempSheet <- pooledSheet
tempSheet$CellLineSuppressed <- factor(substr(tempSheet$CellLine, 3, 4))
tempSheet$Library <- factor(tempSheet$Library, levels=c("Input", "EZH2"), labels=c("Input", "Pulldown"))
tempSheet$Condition <- factor(tempSheet$Condition, levels=c("WT", "MT"), labels=c("WT", "MT"))
```

DESeq2

``` r
tempCount <- pooled.fCounts.count[, tempSheet$Sample]
deseq2.obj <- DESeqDataSetFromMatrix(countData=tempCount, colData=tempSheet, design= ~ Condition + Condition:CellLineSuppressed + Condition:Library)
```

    ## converting counts to integer mode

``` r
sizeFactors(deseq2.obj) <- tempSheet$DeduppedReadNum/min(tempSheet$DeduppedReadNum)
deseq2.obj <- DESeq(deseq2.obj[apply(counts(deseq2.obj),1,max)>=3,])
```

    ## using pre-existing size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
resultsNames(deseq2.obj)
```

    ## [1] "Intercept"                       "Condition_MT_vs_WT"             
    ## [3] "ConditionWT.CellLineSuppressedb" "ConditionMT.CellLineSuppressedb"
    ## [5] "ConditionWT.LibraryPulldown"     "ConditionMT.LibraryPulldown"

``` r
enrich[["Diff"]] <- results(deseq2.obj, contrast=list("ConditionMT.LibraryPulldown", "ConditionWT.LibraryPulldown"))
```

Empirical
p-value

``` r
enrich[["Diff"]] <- assignEmpPvalue(enrich[["Diff"]], main="EZH2: WT vs. MT")
```

    ## Step 1... determine cutoff point
    ## Step 2... estimate parameters of null distribution and eta0
    ## Step 3... compute p-values and estimate empirical PDF/CDF
    ## Step 4... compute q-values and local fdr

![](ChIPseq_EZH2_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

Table

``` r
enrich[["Diff"]] <- deseq2table(enrich[["Diff"]], gencode, 0.1)
```

    ## 
    ## Insig   Sig SigFC 
    ## 37614   189    62

### Combined table

``` r
chipseq.ezh2.table <- gencode

for (i in c("baseMean", "log2FoldChange", "stat", "pvalue", "padj", "sig", "sigFC")) {
  chipseq.ezh2.table[, paste("wt", i, sep=".")] <- enrich[["WT"]][match(chipseq.ezh2.table$gene_id, enrich[["WT"]]$gene_id), i]
  chipseq.ezh2.table[, paste("mt", i, sep=".")] <- enrich[["MT"]][match(chipseq.ezh2.table$gene_id, enrich[["MT"]]$gene_id), i]
  chipseq.ezh2.table[, paste("diff", i, sep=".")] <- enrich[["Diff"]][match(chipseq.ezh2.table$gene_id, enrich[["Diff"]]$gene_id), i]
}
chipseq.ezh2.table <- subset(chipseq.ezh2.table, !is.na(wt.stat) | !is.na(mt.stat) | !is.na(diff.stat))
```

# 4\. Save

``` r
chipseq.ezh2.sampleSheet <- pooledSheet
save(list=c("chipseq.ezh2.sampleSheet", "chipseq.ezh2.table"), file="./ChIPseq/ChIPseq_EZH2.Rdata")
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
