# ChIP-seq pipeline

## 0. Prerequsite

### Softwares
[Samtools](http://www.htslib.org/)

[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[BWA](https://github.com/lh3/bwa)

[Cutadapt](https://cutadapt.readthedocs.io/en/stable/)

[FeatureCounts](http://subread.sourceforge.net/)

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

(1) Run BWA
Step 1
`bwa aln INDEX_prefix -f SAMPLE.sai -t 5 -q 5 -l 32 -k 2 SAMPLE.fastq.trimmed.gz`\
INDEX_prefix: BWA genome index.\
-t N: the number of threads, here assumes using 5 threads.\

Step 2
`bwa samse INDEX_prefix SAMPLE.sai SAMPLE.fastq.trimmed.gz | samtools view -Su - | samtools sort -o SAMPLE.bam -`

(2) BAM index

`samtools index SAMPLE.Aligned.sortedByCoord.out.bam`
