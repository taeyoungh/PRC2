import sys
import subprocess # from subprocess import *
import getopt

'''
1. Description
Macs2 bedgraph has regions out of boundary due to extension of reads.
This program will fix them.

2. Usage
$ chipseq_macs2_bdg_fixer.py -i MACS2_OUTPUT.BDG -g GENOME_SIZE -o OUTPUT_FILE_NAME

(1) Input
GENOME_SIZE must have two columns without header: chrom "\t" length
'''

opts,args = getopt.getopt(sys.argv[1:], "hi:g:o:", ["help", "input=", "genome_size=", "output="]) 

for opt, val in opts:
	if opt=="-h":
		print("Usage: bedGraph_scaler.py -i bedgraph -g genome_size -o output"+"\n")
		continue
	if opt in ("-i", "--input"):
		bdgFileName = val
		continue
	if opt in ("-g", "--genome_size"):
		genomeFileName = val
		continue
	if opt in ("-o", "--output"):
		outputFileName = val
		continue

# ================================
# 1. Read annotation file
# ================================

fin = open(genomeFileName)
lines = fin.readlines()
chromLength = {}
for l in lines:
	l = l.strip().split("\t")
	chromLength[l[0]] = int(l[1])

fin.close()

# ================================
# 2. Fix bedgraph
# ================================

fin = open(bdgFileName)
fout = open(outputFileName, "w")

for l in fin:
	[chrom, start, end, value] = l.strip().split("\t")
	start = int(start)
	end = int(end)

	if start >= chromLength[chrom]:
		continue
	elif end >= chromLength[chrom]:
		end = chromLength[chrom]

	out = [chrom, str(start), str(end), value]
	fout.write("\t".join(out)+"\n")

fin.close()
fout.close()

