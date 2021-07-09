#!/usr/bin/env python3.8

import sys, os
from Bio import SeqIO
from pathlib import Path
import pandas as pd

aceIDtaxIDSpecies = "/datater/wu/data/MSSPE-basic/enterovirus/kmer13/intermediate/aceID-taxID-species.tsv"
genomePath = "/datater/wu/data/test_clustering/all_enterovirus_genomes"

df = pd.read_table(aceIDtaxIDSpecies, sep="\t", index_col=0)
dict_specie_aceID = {k: g["accession"].tolist() for k, g in df.groupby("species")}

for species in dict_specie_aceID:
    spe = species.replace(" ", "_")
    outputDir = f"/datater/wu/data/test_clustering/{spe}/dRep_centroid"
    genomeList = list(map(lambda x: f"{genomePath}/{x}.[0-9].fasta", dict_specie_aceID[species]))
    genomeList = " ".join(genomeList)
    # commande = f"dRep cluster -p 6 -sa 0.95 -nc 0.85 {outputDir} -g {genomeList}"
    commande = f"dRep cluster -p 6 -sa 0.95 -nc 0.85 --clusterAlg centroid {outputDir} -g {genomeList}"
    os.system(commande)
