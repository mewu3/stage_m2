#!/usr/bin/env python3.8
import sys
import re
import os
from collections import defaultdict
from Bio import SeqIO
import glob
from random import sample
import pandas as pd
import dask.dataframe as dd
import numpy as np
import time
#
virus = "enterovirus"
# # minimumLen = 25000 # coronavirus genome length
minimumLen = 7000 # enterovirus genome Length

aceIDtaxID = "/datater/wu/data/ncbiTaxonomy/nucl_gb.accession2taxid"
taxIDLineage = "/datater/wu/data/ncbiTaxonomy/ncbi_lineages_2021-04-23.csv"
datadir = "/datater/wu/data/tmp"
output1 = f"/datater/wu/data/tmp/{virus}.tsv"
output2 = f"/datater/wu/data/tmp/{virus}Curated.fasta"

# dict_virus = {
#     11137: "229E",
#     31631: "OC43",
#     277944: "NL63",
#     290028: "HKU1",
#     694009: "SARS",
#     1335626: "MERS",
# }
# ls_virusAcession = list(dict_virus.keys())

dict_virus = {
    138948: "EA",
    138949: "EB",
    138950: "EC",
    138951: "ED",
    147711: "RA",
    147712: "RB",
    463676: "RC",
}
ls_virusAcession = list(dict_virus.keys())

df_aceIDtaxID = pd.read_csv(aceIDtaxID, sep="\t", header=0, usecols=["accession", "taxid"], dtype={"accession":"object", "taxid":"int64"})


df = df_aceIDtaxID[df_aceIDtaxID["taxid"].isin(ls_virusAcession)]
df["specie"] = df["taxid"].map(dict_virus)
df.to_csv(output1, sep="\t", index=False)
df = pd.read_csv(output1, sep="\t", header=0)

ls_files = []

dict_taxIDaceID = {k: g["accession"].tolist() for k, g in df[["taxid", "accession"]].groupby("taxid")}

for taxID in dict_taxIDaceID:
    specie = dict_virus[taxID]
    specieAceID = f"{datadir}/{specie}.txt"
    specieAceIDOpen = open(f"{datadir}/{specie}.txt", "w")
    for aceID in dict_taxIDaceID[taxID]:
        specieAceIDOpen.write(f"{aceID}\n")
    specieAceIDOpen.close()
    fastaFile = f"{datadir}/{specie}.fasta"
    ls_files.append(fastaFile)
    os.system(f"epost -db nuccore -input {specieAceID} -format acc |efetch -format fasta > {fastaFile}")

min = 99999999999999

ls_newFiles = []

for file in ls_files:
    filtedFileName = os.path.splitext(file)[0] + "_filtered.fasta"
    os.system(f"cat {file}|seqkit seq -m {minimumLen} > {filtedFileName}")
    ls_newFiles.append(filtedFileName)
    seqCount = int(os.popen(f"grep '^>' -c {filtedFileName}").read())
    if seqCount < min:
        min = seqCount

output2Open = open(output2, "a")
for file in ls_newFiles:
    records = [r for r in SeqIO.parse(file, "fasta")]
    sampleRecords = sample(records, min)
    SeqIO.write(sampleRecords, output2Open, "fasta")
output2Open.close()

# try:
#     os.remove(output2)
# except OSError:
#     pass
# #
# file_sarsCov2Acession = "/datater/wu/data/tmp/sarscov2.acc"
# ls_sarsCov2Accession = []
# with open(file_sarsCov2Acession, "r") as f:
#     for line in f:
#         line = line.rstrip("\n")
#         ls_sarsCov2Accession.append(line)
#
# output2Open = open(output2, "a")


# file = "/datater/wu/data/coronavirus.fasta"
# records = [r for r in SeqIO.parse(file, "fasta") if r.id not in ls_sarsCov2Accession]
# SeqIO.write(records, "/datater/wu/data/coronavirus_withoutSarsCov2.fasta", "fasta")
