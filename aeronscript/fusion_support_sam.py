#!/usr/bin/python

import sys
import re

fusionfile = sys.argv[1]
alnsamfile = sys.argv[2]
minlen = int(sys.argv[3])
supportfile = sys.argv[4]
valid_sams_file = sys.argv[5]

cigar_re = re.compile("(\d+)([DHIMNPSX=])")

def aln_length_from_cigar(cigar):
	aln_length = 0
	found = cigar_re.findall(cigar)
	for part in found:
		size = int(part[0])
		operation = part[1]
		if operation in ['D', 'M', 'X', '=']:
			aln_length += size
	return aln_length

fusionpos = {}
support = {}
valid_sams = []

with open(fusionfile) as f:
	while True:
		nameline = f.readline()
		seqline = f.readline()
		if not nameline: break
		name = nameline[1:].strip()
		pos = seqline.index('N')
		fusionpos[name] = pos
		support[name] = 0

with open(alnsamfile) as f:
	for l in f:
		if l[0] == '@': continue
		parts = l.strip().split('\t')
		transcript_name = parts[2]
		if transcript_name in fusionpos:
			aln_start = int(parts[3])
			aln_length = aln_length_from_cigar(parts[5])
			if aln_start > fusionpos[transcript_name] - minlen or aln_start + aln_length < fusionpos[transcript_name] + minlen: continue
			valid_sams.append(l.strip())
			if transcript_name not in support: support[transcript_name] = 0
			support[transcript_name] += 1

pairs = [(fusion, support[fusion]) for fusion in support]
pairs.sort(key=lambda x: -x[1])

with open(supportfile, 'w') as f:
	for fusion in pairs:
		f.write(fusion[0] + "\t" + str(fusion[1]) + "\n")

with open(valid_sams_file, 'w') as f:
	for paf in valid_sams:
		f.write(paf + "\n")
