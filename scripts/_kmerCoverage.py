import os
import sys
import re
from collections import defaultdict
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq

input1 = "/datater/wu/data/enterovirus/enterovirus.uniq"
output1 = "/datater/wu/data/enterovirus/evaluation13/aceID-taxID-species.tsv"

dict_aceIDtaxID=defaultdict(str)

Entrez.email = "wu.meiju@outlook.com"

# acessionID = os.popen(f"cat {input1}|grep '^>'|cut -f1 -d' '|sed 's/>//g'").read()
# taxID = os.popen(f"cat {input1}|grep '^>'|cut -f1 -d' '|sed 's/>//g'|epost -db nuccore|esummary -db nuccore|xtract -pattern DocumentSummary -element Caption,TaxId|cut -f2").read()
# acessionID = acessionID.split("\n")
# taxID = taxID.split("\n")
# combine = [' '.join(x) for x in zip(acessionID, taxID)]
#
# for x in combine:
#     ls = x.split()
#     if len(ls) > 1 :
#         ace = x.split()[0]
#         tax = x.split()[1]
#         dict_aceIDtaxID[ace]=tax

handle = Entrez.efetch(db="taxonomy", id="2152662", mode="text")
records = Entrez.read(handle)
for taxon in records:
    print(taxon)
    if taxon["Rank"] == "no rank":
        for t in taxon["LineageEx"]:
            if t["Rank"] == "species":
                species = t["ScientificName"]
                print(species)
