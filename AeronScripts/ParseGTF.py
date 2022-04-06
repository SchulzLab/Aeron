#1       havana  exon    13221   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002312635"; exon_version "1"; tag "basic"; transcript_support_level "1";

import pandas as pd
from pandas import IndexSlice as IDX
from BCBio import GFF

class ParseGTF():
    def __init__(self, gtf):
        data = []
        f = open(gtf)
        for rec in GFF.parse(f):
            for gene_feature in rec.features:
                for transcript_feature in gene_feature.sub_features:
                    for feature in transcript_feature.sub_features:
                        
                        if feature.type != "exon":
                            continue
                        try:
                            exon_number = feature.qualifiers["exon_number"]
                        except KeyError:
                            exon_number = 1
                        try:
                            data.append(
                                {
                                    "gene id": gene_feature.id, 
                                    "transcript id": transcript_feature.id,
                                    "exon id": feature.id, 
                                    "exon start": feature.location.start, 
                                    "exon end": feature.location.end, 
                                    "chr": rec.id, 
                                    "strand": "+" if feature.strand >= 0 else "-", 
                                    "exon number": exon_number[0]
                                }
                            )
                        except KeyError:
                            continue
        self.data = pd.DataFrame(data).set_index(
            ["gene id", "transcript id", "exon id"], drop=False
        )
        f.close()


    def getExons(self, key):
        try:
            exons = self.data.loc[IDX[:, key, :], :]
        except KeyError:
            try:
                exons = self.data.loc[key, :]
            except:
                exons = 0
        
        return exons
	
    def getStart(self, key):
        return self.data.loc[IDX[:, :, key], "exon start"].iloc[0]

    def getEnd(self, key):
        return self.data.loc[IDX[:, :, key], "exon end"].iloc[0]

    def getChromosome(self, key):
        return self.data.loc[IDX[:, :, key], "chr"].iloc[0]

    def getStrand(self, key):
        return self.data.loc[IDX[:, :, key], "strand"].iloc[0]

    def getExonNumber(self, key):
        return self.data.loc[IDX[:, :, key], "exon number"].iloc[0]

    def getTranscript(self, key):
        return self.data.loc[IDX[:, :, key], "transcript id"]

    def getAllGenes(self):
        return self.data.index.get_level_values("transcript id")

    def getTranscriptPosition(self, key):
        exon_ends = self.data.loc[IDX[:, key, :], "exon end"]
        startpoint, endpoint = exon_ends.min(), exon_ends.max()

        return startpoint,endpoint

    def getAllTranscripts(self):
        return self.data.index.get_level_values("gene id")
