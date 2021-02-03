# 2021/1/3
# Update Source/Sox2/dedupped2.py to handle single end reads.
# Uniquely mapped reads (BWA) are checked for duplicates according to position (1st) + sequence (2nd).

import sys
import string

fin = sys.stdin
fout = sys.stdout
#ferr = sys.stderr

def RevComp(seq):
	complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
	rcseq = seq.translate(complements)[::-1]
	return rcseq

# Iterate alignment records/lines.
umi0 = ""
reads = []
for line in fin:
	if line[0] == "@" : # if header part, skip.
		fout.write(line)
		continue

	line = line.strip().split("\t")

	# Extract umi and umiQual from qname.
	temp = line[0].split(":")
	umiQual = temp[0]
	umi = temp[1]

	# Replace qname with the original qname (without UMI info).
	line[0] = ":".join(temp[3:])[1:] # orginal qname without '@'

	# Skip if UMI has lower sequence quality.
	if umiQual == "FAIL" :
		continue

	# Take uniquely-mapped reads only: BWA
	if line[11] != "XT:A:U":
		continue

	# Start main processing
	if umi0 == "" : # initialization
		umi0 = umi
		reads.append(line)
		continue
	elif umi == umi0 : # collecting the reads with same UMIs
		reads.append(line)
		continue
	else:
		# process the reads with the same UMIs.
		# here, we will classify reads by position (chromosome+location)
		# then, check duplicate of sequence for a given position.

		deduppedIdx = {} # dedupped reads are stored in a dictionary: position is key, (index, seq) is value.
		i = 0
		while i < (len(reads)-1) : # for i-th record with the same umi,

			# find position for the current read.
			chrom = reads[i][2]
			coord = reads[i][3]
			pos = chrom + ":" + coord

			# find sequence for the current read.
			seq = reads[i][9]
			flag = int(reads[i][1])
			if flag & 0b10000 == 0b10000: # reverse strand
				seq = RevComp(seq)

			# checking duplicates
			if pos in deduppedIdx.keys(): # if there is a read with the same position,
				if len([d for d in deduppedIdx[pos] if d[1] != seq]) == len(deduppedIdx[pos]):
					# add only when the current read is new, in other words,
					# add only when all the existing reads are different from the current read.
					deduppedIdx[pos].append([i, seq])
			else: # initialization or add
				deduppedIdx[pos] = [[i, seq]]

			# move to the next read.
			i = i + 1

		# print out the dedupped reads.
		for idxList in deduppedIdx.values():
			for idx in idxList:
				fout.write("\t".join(reads[idx[0]])+"\n")

		# after the current UMI processing is done, update umi.
		umi0 = umi
		reads = [line]

fin.close()
fout.close()
