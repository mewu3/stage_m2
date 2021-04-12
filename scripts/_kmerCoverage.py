import os
import sys
import re
from collections import defaultdict
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq

input1 = "/datater/wu/data/enterovirus/enterovirus.uniq"
input2 = "/datater/wu/data/enterovirus/dimer13/allOligos_reverse.set.fasta"
output1 = "/datater/wu/data/enterovirus/evaluation13/aceID-taxID-species.tsv"

splitFilesDir = "/datater/wu/data/enterovirus/splitFiles"

dict_aceIDtaxID=defaultdict(str)
dict_posiSeq = defaultdict(str)
dict_posiAce = defaultdict(list)
dict_taxIDSpecies = defaultdict(str)
dict_speciesCount = defaultdict(int)

Entrez.email = "wu.meiju@outlook.com"

acessionID = os.popen(f"cat {input1}|grep '^>'|cut -f1 -d' '|sed 's/>//g'").read()
taxID = os.popen(f"cat {input1}|grep '^>'|cut -f1 -d' '|sed 's/>//g'|epost -db nuccore|esummary -db nuccore|xtract -pattern DocumentSummary -element Caption,TaxId|cut -f2").read()
acessionID = acessionID.split("\n")
taxID = taxID.split("\n")
combine = [' '.join(x) for x in zip(acessionID, taxID)]

for x in combine:
    ls = x.split()
    if len(ls) > 1 :
        ace = x.split()[0]
        tax = x.split()[1]
        dict_aceIDtaxID[ace]=tax

output1_open = open(output1, "w")
for aceID in dict_aceIDtaxID:
    taxID = dict_aceIDtaxID[aceID]
    handle = Entrez.efetch(db="taxonomy", id=taxID, mode="text")
    records = Entrez.read(handle)
    for taxon in records:
        for t in taxon["LineageEx"]:
            if t["Rank"] == "species":
                species = t["ScientificName"]
                output1_open.write(f"{aceID}\t{taxID}\t{species}\n")
                if len(species) > 0 :
                    dict_taxIDSpecies[taxID] = species
                else:
                    dict_taxIDSpecies[taxID] = "unknow"
output1_open.close()

# for taxId in dict_taxIDSpecies:
#     species = dict_taxIDSpecies[taxId]
#     dict_speciesCount[species]+=1
#
# input2_open = open(input2, "r")
# for line in input2_open.readlines():
#     line = line.rstrip("\n")
#     if line.startswith(">"):
#         header = line
#         position = header.split("|")[1].split()[1]
#         dict_posiSeq[position] = ""
#     else:
#         kmer = Seq(line)
#         oligo = str(kmer.reverse_complement())
#         dict_posiSeq[position] = oligo
# input2_open.close()
#
# for posi in dict_posiSeq:
#     seq = dict_posiSeq[posi]
#     splitFiles = f"{splitFilesDir}/reverse{posi}.fasta"
#     accessionID_ls = []
#     file = open(splitFiles, "r")
#     for line in file:
#         line = line.rstrip("\n")
#         if line.startswith(">"):
#             accessionID = line.split()[1]
#         if not line.startswith(">"):
#             if re.search(seq, line, re.I):
#                 accessionID_ls.append(accessionID)
#     file.close()
#     dict_posiAce[posi]=accessionID_ls
#
# output1_open = open(output1, "w")
# for posi in dict_posiAce:
#     dict=defaultdict(int)
#     for acessionID in dict_posiAce[posi]:
#         taxID = dict_aceIDtaxID[acessionID]
#         species = dict_taxIDSpecies[taxID]
#         dict[species]+=1
#     for species in dict:
#         species_count = dict[species]
#         totolCount = dict_speciesCount[species]
#         output1_open.write(f"{posi}\t{species}\t{species_count}\t{totolCount}\n")
# output1_open.close()
