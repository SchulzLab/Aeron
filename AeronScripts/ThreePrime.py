from collections import defaultdict
from pprint import pprint
import argparse
import sys

matrixfile = ''
	#opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["ifile=","ofile="])
parser = argparse.ArgumentParser(description='Quantification Counter')
parser.add_argument('-m', action="store", dest="matrix")

par = parser.parse_args()

matrixfile = par.matrix

#Reading the matrix file
m=open(matrixfile)
mt=defaultdict()
mt1 = defaultdict()
max=0
nam=""
tmp="hi"
for line in m:
    ent=line.rstrip().split("\t")
    if(float(ent[-2])>0.2):
        if(ent[0] != nam):
            max=0
            mt[ent[0]]=""
        if(float(ent[-2])>max):
            mt[ent[0]]=[]
            tran=ent[1].split(".")[0]
            mt[ent[0]] = tran
            mt1[ent[0]] = float(ent[-1])
            nam=ent[0]
            max=float(ent[-2])
        elif(float(ent[-1])==max):
            tran=ent[1].split(".")[0]
            if(float(ent[-1])<float(mt1[ent[0]])):
                mt[ent[0]] = tran
                mt1[ent[0]] = float(ent[-1])
        else:
            tmp="Hello"

mtf=defaultdict()
traf = defaultdict()
count=0
for k in mt.keys():
    traf[mt[k]]=0

for reads in mt.keys():
    traf[mt[reads]]+=1

print("Transcript	Count")
for transcript in traf.keys():
    print(transcript+"	"+str(traf[transcript])) 
