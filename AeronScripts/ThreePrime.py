import json
from pprint import pprint
from ParseGTF import *
import argparse
import sys

gtffile = ''
matrixfile = ''
jsonfile = ''
	#opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["ifile=","ofile="])
parser = argparse.ArgumentParser(description='Quantification Counter')
parser.add_argument('-g', action="store", dest="gtf")
parser.add_argument('-m', action="store", dest="matrix")
parser.add_argument('-j', action="store", dest="json")

par = parser.parse_args()

gtffile = par.gtf
matrixfile = par.matrix
jsonfile = par.json

#Reading and analysing gtf file
gf=ParseGTF(gtffile)
readen=defaultdict(list)
readst=defaultdict(list)

#Reading and analysing json file
f= open(jsonfile)
slength=0
teshmap={}
for line in f:
	data = json.loads(line)
	sstart=0
	for i in range (0,len(data['path']['mapping'])):
		node = data['path']['mapping'][i]['position']['name']
		if(sstart==0):
			sstart=(node.split("-"))[1]
		if('sequence' in data['path']['mapping'][i]['edit'][0].keys()):
			slength=len(data['path']['mapping'][i]['edit'][0]['sequence'])
	fnode=(node.split("-"))[1]
	nme = 	data['name']
	readen[nme].append(int(fnode))
	readst[nme].append(int(sstart))

#Reading the matrix file
m=open(matrixfile)
mt=defaultdict(list)
max=0
nam=""
tmp="hi"
for line in m:
	ent=line.rstrip().split("\t")
	if(ent[0] != nam):
		max=0
		mt[ent[0]]=[]
	if(float(ent[-1])>max):
		mt[ent[0]]=[]
		tran=ent[1].split(".")[0]
		mt[ent[0]].append(tran)
		nam=ent[0]
		max=float(ent[-1])
	elif(float(ent[-1])==max):
		tran=ent[1].split(".")[0]
		mt[ent[0]].append(tran)
	else:
		tmp="Hello"	

mtf=defaultdict()
traf = defaultdict()
count=0
for k in mt.keys():
	mt[k] = list(mt[k])
	trns = mt[k][0]
	if(len(mt[k])>1):
		readends = readen[k]
		readstart = readst[k]
		mn=0
		for i in mt[k]:
			transtart, transend = gf.getTranscriptPosition(i)
			if(transtart!=0 and transend!=0):
				for j in range(0,len(readends)):
					dif=transend-readends[j]				
					if(dif>0):
						if(mn==0):
							mn=dif
							trns = i
						if(dif<mn):
							mn = dif
							trns = i
	mtf[k]=trns
	traf[trns]=0

for reads in mtf.keys():
	traf[mtf[reads]]+=1
print("Transcript	Count")
for transcript in traf.keys():
	print(transcript+"	"+str(traf[transcript])) 
