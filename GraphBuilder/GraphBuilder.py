#1       havana  exon    11869   12227   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "1"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002234944"; exon_version "1"; tag "basic"; transcript_support_level "1";
#1       havana  exon    12613   12721   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00003582793"; exon_version "1"; tag "basic"; transcript_support_level "1";
#1       havana  exon    13221   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002312635"; exon_version "1"; tag "basic"; transcript_support_level "1";

import os
import sys
import re
import optparse
from optparse import OptionParser
from optparse import Option, OptionValueError
from ParseGTF import *
        
class ParseOptions():
	def getoptions(self):
		parser = OptionParser()
		wd = os.getcwd()
		parser.add_option('-e', '--genome',dest='ef',help='File containing genome sequences divided into chromosomes', action="store")
		parser.add_option('-g', '--gtf',dest='gf',help='GTF file for the genome', action="store")
		parser.add_option('-o', '--out',dest='out',help='Output Folder', action="store", default=wd)
		(options, args) = parser.parse_args()
		return options

class SequenceAnalyser():
	def FParser(self, seq):
		Sequences = {}
		f=open(seq)
		line = f.readline()
		seq_id=""		
		sequence=""
		while line:
			if (re.match(r'^\>', line.rstrip())):
				Sequences[seq_id] = sequence
				seq_id = (line.rstrip().split(" "))[0]
				sequence=""				
			else:
				sequence=sequence+(line.rstrip())
			line = f.readline()
		Sequences[seq_id] = sequence
		return Sequences
	
	def getExonSequences(self, Sequences, exons, gb):
		Esequences = {}
		for exon in exons:
			key = ">"+gb.getChromosome(exon)[0]
			sequence = Sequences[key]
			length = gb.getEnd(exon) - gb.getStart(exon) 
			Esequence[exon] = sequence[(gb.getStart(exon)-1):gb.getStart(exon)+length]
		return Esequences

class GraphBuild():
	def SortExons(self, gb, exons):
		start = []
		end = []
		chromosome = []
		exn = []
		for exon in exons:
			start.append(gb.getStart(exon))
			chromosome.append(gb.getChromosome(exon)[0])
			end.append(gb.getEnd(exon))
		end = [ed for _, ed in sorted(zip(start,end))]
		chromosome = [cd for _, cd in sorted(zip(start,chromosome))]		
		start = sorted(start)
		return start,end,chromosome
	
	def DefineCluster(self, start, end, chromosome):
		cluster = defaultdict(list)
		cchromosome = defaultdict(list)
		strt = start[0]
		ed = end[0]
		j=0;
		for i in range(0, len(start)):
			if(start[i] >= strt and end[i] <= ed):
				if(start[i] not in cluster[j]):
					cluster[j].append(start[i])
					cchromosome[j].append(chromosome[i])
				if(end[i] not in cluster[j]):
					cluster[j].append(end[i])
					cchromosome[j].append(chromosome[i])
			elif(start[i] >= strt and start[i] <= ed):
				if(start[i] not in cluster[j]):
					cluster[j].append(start[i])
					cchromosome[j].append(chromosome[i])
				if(end[i] not in cluster[j]):
					cluster[j].append(end[i])
					cchromosome[j].append(chromosome[i])
				ed = end[i]
			elif (start[i] > ed):
				j=j+1
				strt = start[i]
				ed = end[i]
				if(start[i] not in cluster[j]):
					cluster[j].append(start[i])
					cchromosome[j].append(chromosome[i])
				if(end[i] not in cluster[j]):
					cluster[j].append(end[i])
					cchromosome[j].append(chromosome[i])
			else:
				print("Something went wrong")
				exit()

		for i in range(0, len(cluster)):
			cluster[i] = sorted(cluster[i])
			cchromosome[i] = sorted(cchromosome[i])
		return cluster,cchromosome

	def getNodeID(self, gb, nodes, nodee, exons):
		nid = {}
		for exon in exons:
			enumber = gb.getExonNumber(exon)
			estart = gb.getStart(exon)
			eend = gb.getEnd(exon)
			key = gb.getTranscript(exon)
			for i in range(0, len(nodes)):
				if(nodes[i]>=estart and nodee[i]<=eend):
					nid[nodes[i]]= key+"-"+str(nodes[i])
		return nid 
		
	def getNodePositions(self, cluster, cchromosome):
		nodes = []
		nodee = []
		chrnod = []
		i=1
		if(len(cluster) == 1 and len(cluster[0]) == 2):
			nodes.append(cluster[0][0])
			nodee.append(cluster[0][1])
			chrnod.append(cchromosome[0][0])
		else:
			for i in range(0, len(cluster)):
				for j in range(0,len(cluster[i])-1):
					nodes.append(cluster[i][j])
					nodee.append(cluster[i][j+1])
					chrnod.append(cchromosome[i][j])
		return nodes,nodee,chrnod
	
	def getNodeConnections(self, nodeid, ndst, nden, chrnod, Sequences, nodes, connections):
		bookkeep = {}
		for i in range(0, len(nodeid)):
			bookkeep[nodeid[ndst[i]]]=0

		for i in range(0, len(nodeid)):
			chromo = chrnod[i]
			key=">"+chromo
			if(key in Sequences.keys()):
				sequence = Sequences[key]			
				length = (nden[i]-ndst[i])+1
				if(len(ndst)>0):
					if(len(ndst)==1):
						sseq = sequence[ndst[i]-1:(ndst[i]+length)]			
					if(i<=(len(ndst)-1)):
						sseq = sequence[ndst[i]-1:(ndst[i]+length)]
					nodes.append("S	"+str(nodeid[ndst[i]])+"	"+sseq)
					bookkeep[str(nodeid[ndst[i]])] = 1
		for i in range(0, len(nodeid)):
			if(bookkeep[str(nodeid[ndst[i]])]==1):
				for j in range(i+1, len(nodeid)):
					if(bookkeep[str(nodeid[ndst[j]])]==1):
						connections.append("L	"+str(nodeid[ndst[i]])+"	+	"+str(nodeid[ndst[j]])+"	+	0M")
		return nodes, connections

	def sanityCheck(self, nodes, connections):
		if len(nodes)==0:
			print("Warning: No sequence information added")
		if len(connections)==0:
			print("Warning: No connection information added")
		

po  = ParseOptions().getoptions()
gb  = GraphBuild()
sq  = SequenceAnalyser() 
fn  = po.gf
seq = po.ef
out = po.out

f = open(out, "w")
nodes = []
connections = []	
		
print("Reading Sequences")
fasta = sq.FParser(seq)
print("Done reading sequences")

print("Reading and processing gtf")
pg = ParseGTF(fn)
print("Done reading gtf")

print("Collecting all the genes")
transcripts = pg.getAllGenes()
print("Building graph for:")

for transcript in transcripts:
	ex = pg.getExons(transcript)
	start,end, chromosome = gb.SortExons(pg, ex)
	cluster, cchromosome = gb.DefineCluster(start, end, chromosome)
	ndst,nden,chrnod = gb.getNodePositions(cluster, cchromosome)	
	nodeid=gb.getNodeID(pg, ndst, nden, ex)
	nodes, connections = gb.getNodeConnections(nodeid, ndst, nden, chrnod, fasta, nodes, connections)

for i in nodes:
	f.write(i+"\n")

for i in connections:
	f.write(i+"\n")

gb.sanityCheck(nodes, connections)
