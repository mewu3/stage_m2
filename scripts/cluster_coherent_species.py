import sys, os
from Bio import SeqIO
from pathlib import Path
import pandas as pd

input = "/datater/wu/data/test_clustering/enterovirus_all.cluster.clstr"
aceIDtaxIDSpecies = "/datater/wu/data/MSSPE-basic/enterovirus/kmer13/intermediate/aceID-taxID-species.tsv"
df = pd.read_table(aceIDtaxIDSpecies, sep="\t", index_col=0)
df = df[["accession", "species"]]
dict = dict(df.values)

dict_cluster = {}

with open(input, "r") as f:
    for l in f:
        l = l.rstrip("\n")
        if l.startswith(">"):
            cluster = l
            dict_cluster[cluster]=set()
        else:
            l=l.split()
            aceID=l[2].lstrip(">").rstrip("...").split(".")[0]
            try:
                specie = dict[aceID]
                dict_cluster[cluster].add(specie)
            except KeyError:
                pass

for cluster in dict_cluster:
    if len(dict_cluster[cluster]) > 1 :
        print(cluster, dict_cluster[cluster])
