#1       havana  exon    11869   12227   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "1"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002234944"; exon_version "1"; tag "basic"; transcript_support_level "1";
#1       havana  exon    12613   12721   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00003582793"; exon_version "1"; tag "basic"; transcript_support_level "1";
#1       havana  exon    13221   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002312635"; exon_version "1"; tag "basic"; transcript_support_level "1";

import os
import re
from tqdm import tqdm
from collections import defaultdict
from optparse import OptionParser
import pandas as pd
from ParseGTF import ParseGTF
        
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
			Esequences[exon] = sequence[(gb.getStart(exon)-1):gb.getStart(exon)+length]
		return Esequences


def GraphBuild(exons: pd.DataFrame, Sequences):
	splice_sites = exons[["exon start", "exon end"]].melt().sort_values("value")
	nodes = getNodePositions(splice_sites)
	nodes = getNodeID(nodes, exons)
	nodes_out, corrections_out = getNodeConnections(nodes, Sequences)
	
	return nodes_out, corrections_out

def getNodePositions(splice_sites: pd.DataFrame):
	rows = splice_sites.iterrows()
	_, prev_site = next(rows)
	nodes = []
	for _, curr_site in rows:
		if prev_site["variable"] == "exon end" and curr_site["variable"] == "exon start":
			nodes.append({"start": prev_site["value"], "end": curr_site["value"]})
		prev_site = curr_site

	return pd.DataFrame(nodes)

def getNodeID(nodes: pd.DataFrame, exons: pd.DataFrame):
	nid = []
	for _, node in nodes.iterrows():
		is_grown = False
		for _, exon in exons.iloc[::-1].iterrows():
			if (node["start"] >= exon["exon start"] or node["end"] <= exon["exon end"]):
				nid.append(f'{exon["transcript id"]}-{node["start"]}')
				is_grown = True
				break
		if not is_grown:
			raise ValueError("Node does not exist in exonic region.")
	nodes.index = nid

	return nodes

def getNodeConnections(nodes: pd.DataFrame, Sequences):
	nodes_out = [
		f'S {idx}	{sequence[node["start"] - 1 : node["end"]]}' 
		for (idx, node), sequence 
		in zip(nodes.iterrows(), Sequences)
	]
	connections_out = []
	nids = nodes.index.to_list()
	for idx, nid_src in enumerate(nids):
		for nid_dest in nids[idx + 1:]:
			connections_out.append(f"L	{nid_src}	+	{nid_dest}	+	0M")

	return nodes_out, connections_out
	

def sanityCheck(nodes, connections):
	if len(nodes)==0:
		print("Warning: No sequence information added")
	if len(connections)==0:
		print("Warning: No connection information added")
		

if __name__ == "__main__":
	po  = ParseOptions().getoptions()
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

	genes = pg.index.unique("gene id")
	print("Building graph")

	for gene in tqdm(genes):
		nodes_of_gene, connections_of_gene = GraphBuild(pg.loc[gene], fasta)
		nodes += nodes_of_gene
		connections += connections_of_gene

	f.write("\n".join(nodes))
	f.write("\n")
	f.write("\n".join(connections))

	sanityCheck(nodes, connections)

	f.close()
