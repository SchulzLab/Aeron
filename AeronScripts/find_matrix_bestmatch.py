#!/usr/bin/python

import sys

matrix_filename = sys.argv[1]
min_goodness_fraction = float(sys.argv[2])
max_fraction_closeness = float(sys.argv[3])

fractions = {}
with open(matrix_filename) as f:
	for l in f:
		parts = l.strip().split('\t')
		readname = parts[0]
		transcriptname = parts[1]
		fraction = float(parts[2])
		if fraction < min_goodness_fraction: continue
		if readname not in fractions: fractions[readname] = []
		fractions[readname].append((transcriptname, fraction))

for read in fractions:
	fractions[read].sort(key = lambda x: -x[1])
	for f in fractions[read]:
		if f[1] >= fractions[read][0][1] * max_fraction_closeness:
			print(read + '\t' + f[0] + '\t' + str(f[1]))
