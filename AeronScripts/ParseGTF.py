#1       havana  exon    13221   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002312635"; exon_version "1"; tag "basic"; transcript_support_level "1";

import re
import pandas as pd
from pandas import IndexSlice as IDX
from BCBio import GFF


attr_regex = re.compile(r'([^\s";]+)[\s=]([^;]+)')


def ParseGTF(gtf):
    data = []
    # with open(gtf, "r") as file:
    #     for row in file:
    #         row.strip().split("\t", maxsplit=8)
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
                    except:
                        continue
    f.close()

    return pd.DataFrame(data).set_index(
        ["gene id", "transcript id", "exon id"], drop=False
    )
