from sklearn import metrics
import sys, os
from Bio import SeqIO
from pathlib import Path
import pandas as pd
import os.path

cdhit = f"/datater/wu/data/test_clustering/all_enterovirus/CD-HIT/enterovirus_all.cluster.clstr"
dRep = f"/datater/wu/data/test_clustering/all_enterovirus/dRep/data_tables/Cdb.csv"
output = "/datater/wu/data/test_clustering/all_enterovirus/compare_CDHIT_dRep.tsv"

dict={}

if os.path.isfile(cdhit):
    with open(cdhit, "r") as f :
        for l in f:
            l = l.rstrip("\n")
            if l.startswith(">"):
                cluster_CDHIT = l.split()[1]
            else:
                genome = l.split()[2].lstrip(">").rstrip("...")
                genome = str(genome.split(".")[0])
                dict[genome] = []
                dict[genome].append(cluster_CDHIT)

if os.path.isfile(dRep):
    with open(dRep, "r") as f:
        for l in f.readlines()[1:]:
            l = l.rstrip("\n").split(",")
            genome = l[0].replace(".fasta", "")
            genome = genome.split(".")[0]
            cluster_dRep = l[1].split("_")[1]
            # print(genome, cluster_dRep)
            dict[genome].append(cluster_dRep)

output = open(output, "w")
output.write(f"Genome\tCD-HIT\tdRep\n")
cluster_CDHIT = []
cluster_dRep = []
for genome in dict:
    cluster_CDHIT.append(dict[genome][0])
    cluster_dRep.append(dict[genome][1])
    output.write(f"{genome}\t{dict[genome][0]}\t{dict[genome][1]}\n")
output.close()

print(metrics.adjusted_rand_score(cluster_CDHIT, cluster_dRep))
print(metrics.adjusted_rand_score(cluster_dRep, cluster_CDHIT))


aceIDtaxIDSpecies = "/datater/wu/data/MSSPE-basic/enterovirus/kmer13/intermediate/aceID-taxID-species.tsv"

df = pd.read_table(aceIDtaxIDSpecies, sep="\t", index_col=0)
dict_specie_aceID = {k: g["accession"].tolist() for k, g in df.groupby("species")}

for species in dict_specie_aceID:
    spe = species.replace(" ", "_")
    dir = f"/datater/wu/data/test_clustering/{spe}"
    cdhit = f"{dir}/CD-HIT/{spe}.cluster.clstr"
    dRep = f"{dir}/dRep/data_tables/Cdb.csv"
    output = f"{dir}/compare_CDHIT_dRep.tsv"

    dict={}

    if os.path.isfile(cdhit):
        with open(cdhit, "r") as f :
            for l in f:
                l = l.rstrip("\n")
                if l.startswith(">"):
                    cluster_CDHIT = l.split()[1]
                else:
                    genome = l.split()[2].lstrip(">").rstrip("...")
                    genome = genome.split(".")[0]
                    dict[genome] = []
                    # print(genome, cluster_CDHIT)
                    dict[genome].append(cluster_CDHIT)

    if os.path.isfile(dRep):
        with open(dRep, "r") as f:
            for l in f.readlines()[1:]:
                l = l.rstrip("\n").split(",")
                genome = l[0].replace(".fasta", "")
                genome = genome.split(".")[0]
                cluster_dRep = l[1].split("_")[1]
                # print(genome, cluster_dRep)
                dict[genome].append(cluster_dRep)

    output = open(output, "w")
    output.write(f"Genome\tCD-HIT\tdRep\n")
    cluster_CDHIT = []
    cluster_dRep = []
    for genome in dict:
        cluster_CDHIT.append(dict[genome][0])
        cluster_dRep.append(dict[genome][1])
        output.write(f"{genome}\t{dict[genome][0]}\t{dict[genome][1]}\n")
    output.close()

    # print(spe)
    # print(metrics.adjusted_rand_score(cluster_CDHIT, cluster_dRep))
    # print(metrics.adjusted_rand_score(cluster_dRep, cluster_CDHIT))
