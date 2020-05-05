# ChIP-seq pipeline

## 0. Prerequsite

### Softwares
[Samtools](http://www.htslib.org/)

[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[BWA](https://github.com/lh3/bwa)

[Cutadapt](https://cutadapt.readthedocs.io/en/stable/)

[PICARD](https://github.com/broadinstitute/picard)

### Genome & annotation
[Gencode vM16](https://www.gencodegenes.org/mouse/release_M16.html)\
Select "comprehensive gene annotation as gtf".

## 1. Fastq

### 1-1. Fastqc
Check quality by running fastqc.

### 1-2. Trim adapter.
`gzip -cd SAMPLE.fastq.gz | cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC - | gzip > SAMPLE.fastq.trimmed.gz`\
SAMPLE.fastq.gz: Sequencing file.\
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC: Illumina adapter.\

## 2. Alignment

### 2-1. BWA

(1) Run BWA\
Step 1\
`bwa aln INDEX_prefix -f SAMPLE.sai -t 5 -q 5 -l 32 -k 2 SAMPLE.fastq.trimmed.gz`\
INDEX_prefix: BWA genome index.\
-t N: the number of threads, here assumes using 5 threads.

Step 2\
`bwa samse INDEX_prefix SAMPLE.sai SAMPLE.fastq.trimmed.gz | samtools view -Su - | samtools sort -o SAMPLE.bam -`

(2) BAM index

`samtools index SAMPLE.bam`

### 2-2. Filtering alignments
Filter unmapped or secondary or failed or duplicated or low-aligned (Phred score<30) or noncanonical chromosomes-aligned reads.

`samtools view -F 1804 -q 30 -b SAMPLE.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM > SAMPLE.filtered.bam`

### 2-3. Removing duplicates
