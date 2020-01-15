#!/usr/bin/python

import sys

matrix_filename = sys.argv[1]
min_valid_fraction = float(sys.argv[2])
max_secondary_fraction = float(sys.argv[3])

count_over_max_secondary = {}
max_fraction = {}
with open(matrix_filename) as f:
	for l in f:
		parts = l.strip().split('\t')
		readname = parts[0]
		transcriptname = parts[1]
		fraction = float(parts[2])
		if fraction < max_secondary_fraction and fraction < min_valid_fraction: continue
		if readname not in max_fraction: max_fraction[readname] = ('', 0)
		if readname not in count_over_max_secondary: count_over_max_secondary[readname] = 0
		if fraction >= max_secondary_fraction: count_over_max_secondary[readname] += 1
		if fraction >= max_fraction[readname][1]: max_fraction[readname] = (transcriptname, fraction)

for read in max_fraction:
	if count_over_max_secondary[read] == 1 and max_fraction[read][1] >= min_valid_fraction:
		print(read + '\t' + max_fraction[read][0] + '\t' + str(max_fraction[read][1]))
