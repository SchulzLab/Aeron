#!/usr/bin/python

import fileinput

pairs = {}

for line in fileinput.input():
	l = line.strip()
	parts = l.split('\t')
	if len(parts) < 3: continue
	readname = '_'.join(parts[0].split('_')[:-2])
	partnum = int(parts[0].split('_')[-2][4:])
	if readname not in pairs: pairs[readname] = {}
	if partnum not in pairs[readname]: pairs[readname][partnum] = ('', '')
	if parts[0][-1] == '1':
		pairs[readname][partnum] = (parts[1], pairs[readname][partnum][1])
	else:
		pairs[readname][partnum] = (pairs[readname][partnum][0], parts[1])

result = {}
current = 0
for r in pairs:
	for n in pairs[r]:
		left = pairs[r][n][0]
		right = pairs[r][n][1]
		if right < left:
			(right, left) = (left, right)
		if (left, right) not in result: result[(left, right)] = {}
		result[(left, right)][r] = True

for p in result:
	for r in result[p]:
		print(p[0] + '\t' + p[1] + '\t' + r)
