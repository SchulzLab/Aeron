#!/usr/bin/python

import fileinput
import re

generegex = re.compile(" gene:(ENSG\d{11}\.\d{1,2}) ")

result = {}

for l in fileinput.input():
	parts = l.strip().split('\t')
	if len(parts) < 3: continue
	lgene = generegex.search(parts[0]).group(1)
	rgene = generegex.search(parts[1]).group(1)
	if lgene == rgene: continue
	if rgene < lgene:
		(rgene, lgene) = (lgene, rgene)
	if (lgene, rgene) not in result: result[(lgene, rgene)] = set()
	result[(lgene, rgene)].add(parts[2])

for p in result:
	print(p[0] + '\t' + p[1] + '\t' + str(len(result[p])))
