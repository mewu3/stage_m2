import os
import sys
import re
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np

oligoSet = "/datater/wu/data/MSSPE-variant/enterovirusCurated/kmer13/allOligo.set"
aceID_taxID_species = "/datater/wu/data/MSSPE-variant/enterovirusCurated/kmer13/intermediate/aceID-taxID-species.tsv"
alignment = "/datater/wu/data/MSSPE-variant/enterovirusCurated/kmer13/intermediate/allKmerCount.sorted.calculated.bowtie"
outputFile = "/datater/wu/data/MSSPE-variant/enterovirusCurated/kmer13/evaluation/genomeCoverage.tsv"

oligoId = []
dict_idPosi = defaultdict(list)
dict_aceIdSpecie = defaultdict(str)
dict_speciesCount = defaultdict(int)

with open(oligoSet, "r") as file:
    for li in file.readlines()[1:]:
        li = li.rstrip("\n").split()
        oligoId.append(f"p{li[0]}")
        dict_idPosi[li[0]] = li[10:]

df_bowtie = pd.read_table(alignment, comment="@", sep="\t", names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "OPT1", "OPT2", "OPT3", "OPT4"])
df_bowtie = df_bowtie[df_bowtie["FLAG"] != 4 ]
df_bowtie = df_bowtie[df_bowtie["QNAME"].isin(oligoId)][["QNAME", "RNAME"]]

dict_idAce = {k: g["RNAME"].tolist() for k, g in df_bowtie.groupby("QNAME")}

for id in dict_idAce:
    dict_idAce[id] = list(map(lambda x: x.split(".")[0], dict_idAce[id]))

input2_open = open(aceID_taxID_species, "r")
for li in input2_open.readlines()[1:]:
    li = li.rstrip("\n").split("\t")
    aceId = li[1]
    specie = li[4]
    dict_aceIdSpecie[aceId]=specie
    dict_speciesCount[specie] += 1
input2_open.close()

dict_sth = defaultdict(int)
for id in dict_idAce:
    for aceID in dict_idAce[id]:
        dict_sth[aceID]+=1

# dict = defaultdict(int)
# for specie in dict_speciesCount:
#     dict[specie]=0

outputOpen = open(outputFile, "w")
for x in range(10):
    dict = defaultdict(int)
    for specie in dict_speciesCount:
        dict[specie]=0
    for ace in dict_sth:
        spe = dict_aceIdSpecie[ace]
        if dict_sth[ace] > x :
            dict[spe]+=1
    for spe in dict:
        species_count = dict[spe]
        totolCount = dict_speciesCount[spe]
        outputOpen.write(f"{x+1}\t{spe}\t{species_count}\t{totolCount}\n")
outputOpen.close()

# for ace in dict_sth:
#     spe = dict_aceIdSpecie[ace]
#     for x in range(10):
#         if dict_sth[ace] > x :
#             dict[spe]+=1
#     for spe in dict:
#         species_count = dict[spe]
#         totolCount = dict_speciesCount[spe]
#         outputOpen.write(f"{x}\t{spe}\t{species_count}\t{totolCount}\n")
