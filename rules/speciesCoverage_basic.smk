#!/usr/bin/env python3.8
import os

configfile: "config.yaml"

samples = list(config["samples"].keys())
kmerSize = list(config["kmerSize"].keys())

file_refSeq = config["file_refSeq"]
file_aceIDtaxID = config["file_aceIDtaxID"]
file_taxIDLineage = config["file_taxIDLineage"]

dataDir = config["dataDir"]

deduplication = config["deduplication"]
clustering = config["clustering"]

clusterIdentity = config["cd-hit"]["identity"]

rule all: 
    input: 
        expand(
            f"{dataDir}/{{sample}}/{{kmerSize}}/evaluation/allOligo.set.totalCoverage5",
            sample = samples,
            kmerSize = kmerSize
        )

rule evaluation2:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/allOligo.set",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/aceID-taxID-species.tsv"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/evaluation/allOligo.set.totalCoverage5"
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
            # print(posi, accessionID_ls)

        input2_open = open(input2, "r")
        for li in input2_open.readlines()[1:]:
            li = li.rstrip("\n").split("\t")
            aceId = li[1]
            specie = li[4]
            dict_aceIdSpecie[aceId]=specie
            dict_speciesCount[specie] += 1
        input2_open.close()

        dict_sth = defaultdict(int)
        for posi in dict_posiAce:
            for aceID in dict_posiAce[posi]: 
                dict_sth[aceID] += 1 

        dict = defaultdict(int) 
        for specie in dict_speciesCount:
            dict[specie]=0
        
        for ace in dict_sth:  
            spe = dict_aceIdSpecie[ace]   
            if dict_sth[ace] > 5 : 
                dict[spe]+=1
        print(dict)
        
        output1_open = open(output1, "w")
        for spe in dict: 
            species_count = dict[spe]
            totolCount = dict_speciesCount[spe]
            output1_open.write(f"{spe}\t{species_count}\t{totolCount}\n")
        output1_open.close()
        
