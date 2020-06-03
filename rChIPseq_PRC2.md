rChIP-seq analysis - PRC2 (EZH2 and SUZ12) enrichment
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
sampleSheet <- subset(sampleSheet, Experiment =="rChIPseq_PRC2" | (Experiment=="ChIPseq" & Library=="Lysate_Input" & CellLine %in% c("WT_A", "WT_B") & Treatment=="None"))
```

### featureCounts output

``` r
ip.fCounts.table <- read.table("./ChIPseq/rChIPseqPRC2_IP_fCounts.txt", header=T, sep="\t", stringsAsFactors = F, check.names = F)
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

Prepare pooling
table

``` r
sampleSheet$CellLine <- factor(sampleSheet$CellLine, levels=c("WT_A", "WT_B"), labels=c("WTa", "WTb"))
sampleSheet$Treatment <- factor(sampleSheet$Treatment, levels=c("RNaseA_neg", "RNaseA_pos", "None"), labels=c("RNaseAneg", "RNaseApos", "None"))
sampleSheet$Pulldown <- sapply(strsplit(sampleSheet$Library, split="_"), "[[", 2)
sampleSheet$PooledName <- paste(sampleSheet$CellLine, sampleSheet$Treatment, sampleSheet$Pulldown, sep="_")
```

Perform
pooling

``` r
pooled.fCounts.count <- sapply(split(sampleSheet, sampleSheet$PooledName), function(x) return(rowSums(fCounts.count[, x$Sample, drop=F])))
```

### Pooled sheet

``` r
pooledSheet <- data.frame(CellLine=rep(c("WTa", "WTb"), each=4), Treatment=rep(c("RNaseAneg", "RNaseApos"), 4), Pulldown=rep(c("EZH2", "EZH2", "SUZ12", "SUZ12"), 2), stringsAsFactors = F)
pooledSheet$PooledName <- paste(pooledSheet$CellLine, pooledSheet$Treatment, pooledSheet$Pulldown, sep="_")
pooledSheet <- rbind(pooledSheet[, c("PooledName", "CellLine", "Treatment", "Pulldown")], subset(sampleSheet, Pulldown=="Input")[, c("PooledName", "CellLine", "Treatment", "Pulldown")])
```

Add dedupped read number for
normalization

``` r
temp <- dcast(sampleSheet, PooledName~., value.var = "DeduppedNum", fun.aggregate = sum)
pooledSheet$DeduppedReadNum <- temp[match(pooledSheet$PooledName, temp$PooledName), 2]
```

``` r
colnames(pooledSheet) <- c("Sample", "CellLine", "Treatment", "Library", "DeduppedReadNum")
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

# 3-1. EZH2 enrichment

``` r
ezh2.enrich<-list()
```

### RNaseA -

Define
factors

``` r
tempSheet <- subset(pooledSheet, Library %in% c("EZH2", "Input") & Treatment %in% c("RNaseAneg", "None"))
tempSheet$CellLine <- factor(tempSheet$CellLine)
tempSheet$Library <- factor(tempSheet$Library, levels=c("Input", "EZH2"), labels=c("Input", "Pulldown"))
```

DEseq2

``` r
ezh2.enrich[["RNaseAneg"]] <- runDESeq2ForChIPseq(tempSheet, pooled.fCounts.count, " ~ CellLine+Library", normVector=tempSheet$DeduppedReadNum)
```

    ## converting counts to integer mode

    ## using pre-existing size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

Empirical
p-value

``` r
ezh2.enrich[["RNaseAneg"]] <- assignEmpPvalue(ezh2.enrich[["RNaseAneg"]], main="EZH2: RNase A-")
```

    ## Step 1... determine cutoff point
    ## Step 2... estimate parameters of null distribution and eta0
    ## Step 3... compute p-values and estimate empirical PDF/CDF
    ## Step 4... compute q-values and local fdr

![](rChIPseq_PRC2_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

Table

``` r
ezh2.enrich[["RNaseAneg"]] <- deseq2table(ezh2.enrich[["RNaseAneg"]], gencode, 0.05)
```

    ## 
    ## Insig   Sig SigFC 
    ## 34153  2417   145

### RNaseA +

Define
factors

``` r
tempSheet <- subset(pooledSheet, Library %in% c("EZH2", "Input") & Treatment %in% c("RNaseApos", "None"))
tempSheet$CellLine <- factor(tempSheet$CellLine)
tempSheet$Library <- factor(tempSheet$Library, levels=c("Input", "EZH2"), labels=c("Input", "Pulldown"))
```

DEseq2

``` r
ezh2.enrich[["RNaseApos"]] <- runDESeq2ForChIPseq(tempSheet, pooled.fCounts.count, " ~ CellLine+Library", normVector=tempSheet$DeduppedReadNum)
```

    ## converting counts to integer mode

    ## using pre-existing size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

Empirical
p-value

``` r
ezh2.enrich[["RNaseApos"]] <- assignEmpPvalue(ezh2.enrich[["RNaseApos"]], main="EZH2: RNase A+")
```

    ## Step 1... determine cutoff point
    ## Step 2... estimate parameters of null distribution and eta0
    ## Step 3... compute p-values and estimate empirical PDF/CDF
    ## Step 4... compute q-values and local fdr

![](rChIPseq_PRC2_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

Table

``` r
ezh2.enrich[["RNaseApos"]] <- deseq2table(ezh2.enrich[["RNaseApos"]], gencode, 0.05)
```

    ## 
    ## Insig   Sig 
    ## 36668   134

### Differential enrichment

Define factors

``` r
tempSheet <- subset(pooledSheet, Library =="EZH2")
tempSheet$CellLine <- factor(tempSheet$CellLine)
tempSheet$Treatment <- factor(tempSheet$Treatment)
```

DESeq2

``` r
tempCount <- pooled.fCounts.count[, tempSheet$Sample]

deseq2.obj <- DESeqDataSetFromMatrix(countData=tempCount, colData=tempSheet, design=~CellLine+Treatment)
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

    ## [1] "Intercept"          "CellLineWTa"        "CellLineWTb"       
    ## [4] "TreatmentRNaseAneg" "TreatmentRNaseApos"

``` r
ezh2.enrich[["Diff"]] <- results(deseq2.obj, contrast=list("TreatmentRNaseApos", "TreatmentRNaseAneg"))
```

Empirical
p-value

``` r
ezh2.enrich[["Diff"]] <- assignEmpPvalue(ezh2.enrich[["Diff"]], main="EZH2: RNase A- vs. RNase A+")
```

    ## Step 1... determine cutoff point
    ## Step 2... estimate parameters of null distribution and eta0
    ## Step 3... compute p-values and estimate empirical PDF/CDF
    ## Step 4... compute q-values and local fdr

![](rChIPseq_PRC2_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

Table

``` r
ezh2.enrich[["Diff"]] <- deseq2table(ezh2.enrich[["Diff"]], gencode, 0.1)
```

    ## 
    ## Insig   Sig 
    ## 35512   509

### Combined table

``` r
ezh2.table <- gencode

for (i in c("baseMean", "log2FoldChange", "stat", "pvalue", "padj", "sig", "sigFC")) {
  ezh2.table[, paste("rnaseAneg", i, sep=".")] <- ezh2.enrich[["RNaseAneg"]][match(ezh2.table$gene_id, ezh2.enrich[["RNaseAneg"]]$gene_id), i]
  ezh2.table[, paste("rnaseApos", i, sep=".")] <- ezh2.enrich[["RNaseApos"]][match(ezh2.table$gene_id, ezh2.enrich[["RNaseApos"]]$gene_id), i]
  ezh2.table[, paste("diff", i, sep=".")] <- ezh2.enrich[["Diff"]][match(ezh2.table$gene_id, ezh2.enrich[["Diff"]]$gene_id), i]
}
ezh2.table <- subset(ezh2.table, !is.na(rnaseAneg.stat) | !is.na(rnaseApos.stat) | !is.na(diff.stat))
dim(ezh2.table)
```

    ## [1] 36982    28

# 3-2. SUZ12 enrichment

``` r
suz12.enrich<-list()
```

### RNaseA -

Define
factors

``` r
tempSheet <- subset(pooledSheet, Library %in% c("SUZ12", "Input") & Treatment %in% c("RNaseAneg", "None"))
tempSheet$CellLine <- factor(tempSheet$CellLine)
tempSheet$Library <- factor(tempSheet$Library, levels=c("Input", "SUZ12"), labels=c("Input", "Pulldown"))
```

DEseq2

``` r
suz12.enrich[["RNaseAneg"]] <- runDESeq2ForChIPseq(tempSheet, pooled.fCounts.count, " ~ CellLine+Library", normVector=tempSheet$DeduppedReadNum)
```

    ## converting counts to integer mode

    ## using pre-existing size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

Empirical
p-value

``` r
suz12.enrich[["RNaseAneg"]] <- assignEmpPvalue(suz12.enrich[["RNaseAneg"]], main="SUZ12: RNase A-")
```

    ## Step 1... determine cutoff point
    ## Step 2... estimate parameters of null distribution and eta0
    ## Step 3... compute p-values and estimate empirical PDF/CDF
    ## Step 4... compute q-values and local fdr

![](rChIPseq_PRC2_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

Table

``` r
suz12.enrich[["RNaseAneg"]] <- deseq2table(suz12.enrich[["RNaseAneg"]], gencode, 0.05)
```

    ## 
    ## Insig   Sig SigFC 
    ## 35934   554   182

### RNaseA +

Define
factors

``` r
tempSheet <- subset(pooledSheet, Library %in% c("SUZ12", "Input") & Treatment %in% c("RNaseApos", "None"))
tempSheet$CellLine <- factor(tempSheet$CellLine)
tempSheet$Library <- factor(tempSheet$Library, levels=c("Input", "SUZ12"), labels=c("Input", "Pulldown"))
```

DEseq2

``` r
suz12.enrich[["RNaseApos"]] <- runDESeq2ForChIPseq(tempSheet, pooled.fCounts.count, " ~ CellLine+Library", normVector=tempSheet$DeduppedReadNum)
```

    ## converting counts to integer mode

    ## using pre-existing size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

Empirical
p-value

``` r
suz12.enrich[["RNaseApos"]] <- assignEmpPvalue(suz12.enrich[["RNaseApos"]], main="SUZ12: RNase A+")
```

    ## Step 1... determine cutoff point
    ## Step 2... estimate parameters of null distribution and eta0
    ## Step 3... compute p-values and estimate empirical PDF/CDF
    ## Step 4... compute q-values and local fdr

![](rChIPseq_PRC2_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

Table

``` r
suz12.enrich[["RNaseApos"]] <- deseq2table(suz12.enrich[["RNaseApos"]], gencode, 0.05)
```

    ## 
    ## Insig   Sig 
    ## 35964   614

### Differential enrichment

Define factors

``` r
tempSheet <- subset(pooledSheet, Library =="SUZ12")
tempSheet$CellLine <- factor(tempSheet$CellLine)
tempSheet$Treatment <- factor(tempSheet$Treatment)
```

DESeq2

``` r
tempCount <- pooled.fCounts.count[, tempSheet$Sample]

deseq2.obj <- DESeqDataSetFromMatrix(countData=tempCount, colData=tempSheet, design=~CellLine+Treatment)
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

    ## [1] "Intercept"          "CellLineWTa"        "CellLineWTb"       
    ## [4] "TreatmentRNaseAneg" "TreatmentRNaseApos"

``` r
suz12.enrich[["Diff"]] <- results(deseq2.obj, contrast=list("TreatmentRNaseApos", "TreatmentRNaseAneg"))
```

Empirical
p-value

``` r
suz12.enrich[["Diff"]] <- assignEmpPvalue(suz12.enrich[["Diff"]], main="SUZ12: RNase A- vs. RNase A+")
```

    ## Step 1... determine cutoff point
    ## Step 2... estimate parameters of null distribution and eta0
    ## Step 3... compute p-values and estimate empirical PDF/CDF
    ## Step 4... compute q-values and local fdr

![](rChIPseq_PRC2_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

Table

``` r
suz12.enrich[["Diff"]] <- deseq2table(suz12.enrich[["Diff"]], gencode, 0.1)
```

    ## 
    ## Insig   Sig SigFC 
    ## 34271   809     1

### Combined table

``` r
suz12.table <- gencode

for (i in c("baseMean", "log2FoldChange", "stat", "pvalue", "padj", "sig", "sigFC")) {
  suz12.table[, paste("rnaseAneg", i, sep=".")] <- suz12.enrich[["RNaseAneg"]][match(suz12.table$gene_id, suz12.enrich[["RNaseAneg"]]$gene_id), i]
  suz12.table[, paste("rnaseApos", i, sep=".")] <- suz12.enrich[["RNaseApos"]][match(suz12.table$gene_id, suz12.enrich[["RNaseApos"]]$gene_id), i]
  suz12.table[, paste("diff", i, sep=".")] <- suz12.enrich[["Diff"]][match(suz12.table$gene_id, suz12.enrich[["Diff"]]$gene_id), i]
}
suz12.table <- subset(suz12.table, !is.na(rnaseAneg.stat) | !is.na(rnaseApos.stat) | !is.na(diff.stat))
dim(suz12.table)
```

    ## [1] 36803    28

# 4\. Save

``` r
rchipseqPRC2.ezh2.table <- ezh2.table
rchipseqPRC2.suz12.table <- suz12.table
rchipseqPRC2.sampleSheet <- pooledSheet
save(list=c("rchipseqPRC2.ezh2.table", "rchipseqPRC2.suz12.table", "rchipseqPRC2.sampleSheet"), file="./ChIPseq/rChIPseq_PRC2.Rdata")
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
