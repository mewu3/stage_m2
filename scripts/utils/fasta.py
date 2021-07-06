#!/usr/bin/env python3.8

import sys, os
from Bio import SeqIO
from pathlib import Path
import pandas as pd

def splitMultiFasta(multiFasta, outputDir):
    Path(outputDir).mkdir(parents=True, exist_ok=True)
    with open(multiFasta) as multiFasta:
        records = SeqIO.parse(multiFasta, "fasta")
        for record in records:
            outputFile = f"{outputDir}/{record.id}.fasta"
            with open(outputFile, "w") as outputFile:
                SeqIO.write(record, outputFile, "fasta")


def seperateSpecies(multiFasta, aceIDtaxIDSpecies, outputDir):

    df = pd.read_table(aceIDtaxIDSpecies, sep="\t", index_col=0)
    dict_specie_aceID = {k: g["accession"].tolist() for k, g in df.groupby("species")}

    dict_seq = {}

    with open(multiFasta, "r") as multiFasta :
        for l in multiFasta :
            l=l.lstrip("\n")
            if l.startswith(">"):
                header = l
                aceID = header.split("|")[0].split(".")[0].lstrip(">")
                dict_seq[aceID]={}
                dict_seq[aceID][header]=""
            else :
                dict_seq[aceID][header]+=l

    for specie in dict_specie_aceID:
        Path(f"{outputDir}/{specie}").mkdir(parents=True, exist_ok=True)
        for aceID in dict_specie_aceID[specie]:
            for key, value in dict_seq[aceID].items():
                with open(f"{outputDir}/{specie}/{aceID}.fasta", "w") as outputFile:
                    outputFile.write(key)
                    outputFile.write(value)

# seperateSpecies(sys.argv[1], sys.argv[2], sys.argv[3])
