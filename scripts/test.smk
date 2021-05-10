#!/usr/bin/env python3.8
import os

configfile: "config.yaml"
samples = list(config["samples"].keys())
file_refSeq = config["file_refSeq"]
file_aceIDtaxID = config["file_aceIDtaxID"]
file_taxIDLineage = config["file_taxIDLineage"]
dataDir = config["dataDir"]
kmerSize = list(config["kmerSize"].keys())
clusterIdentity = config["cd-hit"]["identity"]

rule all:
    input:
        expand(
            f"{dataDir}/{{sample}}/{{kmerSize}}/evaluation/allOligo.set.coverage",
            sample = samples,
            kmerSize = kmerSize
        )

rule evaluation2:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/allOligo.set",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/aceID-taxID-species.tsv",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.bowtie"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/evaluation/allOligo.set.coverage"
    run:
        import os
        import sys
        import re
        from collections import defaultdict
        from Bio import SeqIO
        from Bio.Seq import Seq
        import pandas as pd
        import numpy as np

        oligoId = []
        dict_idPosi = defaultdict(list)
        # dict_idSpecieCount = defaultdict(lambda: defaultdict(int))
        dict_aceIdSpecie = defaultdict(str)
        dict_speciesCount = defaultdict(int)

        with open(input[0], "r") as file:
            for li in file.readlines()[1:]:
                li = li.rstrip("\n").split()
                oligoId.append(f"p{li[0]}")
                # position = f"{li[-2]}-{li[-1]}"
                dict_idPosi[li[0]] = li[10:]

        df_bowtie = pd.read_table(input[2], comment="@", sep="\t", names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "OPT1", "OPT2", "OPT3", "OPT4"])
        df_bowtie
        df_bowtie = df_bowtie[df_bowtie["QNAME"].isin(oligoId)][["QNAME", "RNAME"]]

        dict_idAce = {k: g["RNAME"].tolist() for k, g in df_bowtie.groupby("QNAME")}

        for id in dict_idAce:
            dict_idAce[id] = list(map(lambda x: x.split(".")[0], dict_idAce[id]))

        input2_open = open(input[1], "r")
        for li in input2_open.readlines()[1:]:
            li = li.rstrip("\n").split("\t")
            aceId = li[1]
            specie = li[4]
            dict_aceIdSpecie[aceId]=specie
            dict_speciesCount[specie] += 1
        input2_open.close()

        outputOpen = open(output[0], "w")
        for id in dict_idAce:
            dict = defaultdict(int)
            for specie in dict_speciesCount:
                dict[specie]=0
            for ace in dict_idAce[id]:
                specie = dict_aceIdSpecie[ace]
                dict[specie]+=1
            id = id.lstrip("p")
            posi = dict_idPosi[id]
            posi = "\t".join(map(str, dict_idPosi[id]))
            for specie in dict:
                species_count = dict[specie]
                totolCount = dict_speciesCount[specie]
                outputOpen.write(f"{posi}\t{specie}\t{species_count}\t{totolCount}\n")
        outputOpen.close()
