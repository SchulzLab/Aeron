#!/usr/bin/python

import sys

fusion_file = sys.argv[1]
corrected_read_file = sys.argv[2]

exemplars = {}
support = {}

with open(fusion_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if len(parts) < 1: continue
		key = (parts[6], parts[7], parts[8], parts[9])
		lenvalue = min(int(parts[5]), int(parts[10]))
		readname = parts[0]
		if key not in support: support[key] = 0
		support[key] += 1
		if key not in exemplars or lenvalue > exemplars[key][0]:
			exemplars[key] = (lenvalue, readname, (parts[3], parts[5], parts[4], parts[10]))

picked_reads = {exemplars[v][1]: (exemplars[v][2][0], exemplars[v][2][1], exemplars[v][2][2], exemplars[v][2][3], str(support[v])) for v in exemplars}
fusion_num = 1

with open(corrected_read_file) as f:
	while True:
		nameline = f.readline()
		seq = f.readline()
		if not nameline: break
		name = nameline[1:].strip()
		if name in picked_reads:
			print(">fusion_" + str(fusion_num) + "_" + picked_reads[name][0] + "_" + picked_reads[name][1] + "bp_" + picked_reads[name][2] + "_" + picked_reads[name][3]+ "bp_" + picked_reads[name][4] + "reads")
			print(seq.strip())
			fusion_num += 1
