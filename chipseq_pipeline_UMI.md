# ChIP-seq pipeline with UMI

Modify ChIP-seq pipeline to use UMI sequencing.

## 1. Fastq

### 1-1. Fastqc
Check quality by running fastqc.

### 1-2. Put UMI information to Read 1 file
`umiToQname3 Read1.fastq.gz Read2.fastq.gz SAMPLE`\
This generates SAMPLE.fastq.gz which has UMI in header.

### 1-3. Trim adapter.
`gzip -cd SAMPLE.fastq.gz | cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC - | gzip > SAMPLE.fastq.trimmed.gz`\
SAMPLE.fastq.gz: Sequencing file.\
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC: Illumina adapter.

## 2. Alignment

### 2-1. BWA

#### (1) BWA: step1
`bwa aln INDEX_prefix -f SAMPLE.sai -t 5 -q 5 -l 32 -k 2 SAMPLE.fastq.trimmed.gz`\
INDEX_prefix: BWA genome index.\
-t N: the number of threads, here assumes using 5 threads.

#### (2) BWA: step 2
`bwa samse INDEX_prefix SAMPLE.sai SAMPLE.fastq.trimmed.gz | samtools view -Su - | samtools sort -o SAMPLE.bam -`

#### (3) BAM index
`samtools index SAMPLE.bam`

### 2-2. Filtering alignments
Filter unmapped or secondary or failed or duplicated or low-aligned (Phred score<30) or noncanonical chromosomes-aligned reads.

`samtools view -F 1804 -q 30 -b SAMPLE.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM > SAMPLE.filtered.bam`

### 2-3. Removing duplicates by using UMI
First, sort BAM by name.

`samtools sort -n -@5 SAMPLE.filtered.bam -o SAMPLE.nameSorted.bam`

Second, remove duplicated reads. Note that this only handles uniquely-mapped reads.

`samtools view -h SAMPLE.nameSorted.bam | python dedupped_single.py | samtools view -Sb - | samtools sort -@5 - > SAMPLE.uniq.dedupped.bam`

## 3. Peak calling
Follow standard ChIP-seq pipeline.
