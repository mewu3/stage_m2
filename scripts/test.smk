# !/usr/bin/env python3.8
import os

configfile: "config.yaml"
samples = list(config["samples"].keys())
file_refSeq = config["file_refSeq"]
file_aceIDtaxID = config["file_aceIDtaxID"]
file_taxIDLineage = config["file_taxIDLineage"]
dataDir = config["dataDir"]
kmerSize = config["kmerSize"]


## msspe basic
rule all:
    input:
        expand(
            f"{dataDir}/{{sample}}/kmer{kmerSize}/evaluation/allKmerCount.sorted.calculated.filtered.allOligo.set.coverage",
            sample = samples
        )

rule basic:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/evaluation/allOligo_after2.tsv",
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/aceID-taxID-species.tsv"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/evaluation/allKmerCount.sorted.calculated.filtered.allOligo.set.coverage"
    params:
        splitFilesDir = f"{dataDir}/{{sample}}/splitFiles"
    run:
        import os
        import sys
        import re
        from collections import defaultdict
        from Bio import SeqIO
        from Bio.Seq import Seq

        input1 = input[0]
        input2 = input[1]
        output1 = output[0]

        dict_posiSeq = defaultdict(str)
        dict_posiAce = defaultdict(list)
        dict_aceIdSpecie = defaultdict(str)
        dict_speciesCount = defaultdict(int)

        input1_open = open(input1, "r")
        for line in input1_open.readlines()[1:]:
            line = line.rstrip("\n").split()
            position = f"{line[-2]}-{line[-1]}"
            oligo = line[1]
            dict_posiSeq[position] = oligo
        input1_open.close()

        for posi in dict_posiSeq:
            seq = str(Seq(dict_posiSeq[posi]).reverse_complement())
            splitFiles = f"{params.splitFilesDir}/reverse{posi}.fasta"
            accessionID_ls = []
            record_dict = SeqIO.to_dict(SeqIO.parse(splitFiles, "fasta"))
            for record in record_dict:
                if re.search(seq, str(record_dict[record].seq), re.I):
                    acessionID = record_dict[record].id.split(".")[0]
                    accessionID_ls.append(acessionID)
            dict_posiAce[posi]=accessionID_ls

        input2_open = open(input2, "r")
        for li in input2_open.readlines()[1:]:
            li = li.rstrip("\n").split("\t")
            aceId = li[1]
            specie = li[4]
            dict_aceIdSpecie[aceId]=specie
            dict_speciesCount[specie] += 1
        input2_open.close()

        output1_open = open(output1, "w")
        for posi in dict_posiAce:
            dict=defaultdict(int)
            for specie in dict_speciesCount:
                dict[specie]=0
            for acessionID in dict_posiAce[posi]:
                spec = dict_aceIdSpecie[acessionID]
                dict[spec]+=1
            for specie in dict:
                species_count = dict[specie]
                totolCount = dict_speciesCount[specie]
                posi = "\t".join(map(str, posi.split("-")))
                output1_open.write(f"{posi}\t{specie}\t{species_count}\t{totolCount}\n")
        output1_open.close()

# rule variant:
#     input:
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/evaluation/allOligo_after2.tsv",
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/aceID-taxID-species.tsv",
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.bowtie"
#     output:
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/evaluation/allKmerCount.sorted.calculated.filtered.allOligo.set.coverage"
#     run:
#         import os
#         import sys
#         import re
#         from collections import defaultdict
#         from Bio import SeqIO
#         from Bio.Seq import Seq
#         import pandas as pd
#
#         oligoId = []
#         dict_idPosi = defaultdict(list)
#         # dict_idSpecieCount = defaultdict(lambda: defaultdict(int))
#         dict_aceIdSpecie = defaultdict(str)
#         dict_speciesCount = defaultdict(int)
#
#         with open(input[0], "r") as file:
#             for li in file.readlines()[1:]:
#                 li = li.rstrip("\n").split()
#                 oligoId.append(f"p{li[0]}")
#                 # position = f"{li[-2]}-{li[-1]}"
#                 dict_idPosi[li[0]] = li[-4:]
#
#         df_bowtie = pd.read_table(input[2], comment="@", sep="\t", names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "OPT1", "OPT2", "OPT3", "OPT4"])
#
#         df_bowtie = df_bowtie[df_bowtie["QNAME"].isin(oligoId)][["QNAME", "RNAME"]]
#
#         dict_idAce = {k: g["RNAME"].tolist() for k, g in df_bowtie.groupby("QNAME")}
#
#         for id in dict_idAce:
#             dict_idAce[id] = list(map(lambda x: x.split(".")[0], dict_idAce[id]))
#
#         input2_open = open(input[1], "r")
#         for li in input2_open.readlines()[1:]:
#             li = li.rstrip("\n").split("\t")
#             aceId = li[1]
#             specie = li[4]
#             dict_aceIdSpecie[aceId]=specie
#             dict_speciesCount[specie] += 1
#         input2_open.close()
#
#         outputOpen = open(output[0], "w")
#         for id in dict_idAce:
#             dict = defaultdict(int)
#             for specie in dict_speciesCount:
#                 dict[specie]=0
#             for ace in dict_idAce[id]:
#                 specie = dict_aceIdSpecie[ace]
#                 dict[specie]+=1
#             id = id.lstrip("p")
#             posi = dict_idPosi[id]
#             posi = "\t".join(map(str, dict_idPosi[id]))
#             for specie in dict:
#                 species_count = dict[specie]
#                 totolCount = dict_speciesCount[specie]
#                 outputOpen.write(f"{posi}\t{specie}\t{species_count}\t{totolCount}\n")
#         outputOpen.close()
