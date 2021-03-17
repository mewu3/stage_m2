#!/usr/bin/env python3.8
import os

configfile: "config.yaml"

samples = list(config["samples"].keys())
dataDir = config["dataDir"]

def fetchConfigParameters():
    if config["mafft"]["fmodel"] == "on":
        config["mafft"]["fmodel"] = "--fmodel"
    else:
        config["mafft"]["fmodel"] = ""

    if config["mafft"]["clustalout"] == "on":
        config["mafft"]["clustalout"] = "--clustalout"
    else:
        config["mafft"]["clustalout"] = ""

    if config["mafft"]["inputorder"] == "on":
        config["mafft"]["inputorder"] = "--inputorder"
    else :
        config["mafft"]["inputorder"] = ""

    if config["mafft"]["reorder"] == "on":
        config["mafft"]["reorder"] = "--reorder"
    else:
        config["mafft"]["reorder"] = ""

    if config["mafft"]["treeout"] == "on":
        config["mafft"]["treeout"] = "--treeout"
    else :
        config["mafft"]["treeout"] = ""

    if config["mafft"]["quiet"] == "on":
        config["mafft"]["quiet"] = "--quiet"
    else:
        config["mafft"]["quiet"] = ""
fetchConfigParameters()

# def aggregate_forwardInput(wildcards):
#     checkpoint_output = checkpoints.splitFiles.get(**wildcards).output[0]
#     return expand("{}/{{sample}}/splitFiles/forward{{seg}}.fasta".format(dataDir),
#                   sample = samples,
#                   seg = glob_wildcards(os.path.join(checkpoint_output, "forward{seg}.fasta")).seg)
#
# def aggregate_reverseInput(wildcards):
#     checkpoint_output = checkpoints.splitFiles.get(**wildcards).output[0]
#     return expand("{}/{{sample}}/splitFiles/reverse{{seg}}.fasta".format(dataDir),
#                   sample = samples,
#                   seg = glob_wildcards(os.path.join(checkpoint_output, "reverse{seg}.fasta")).seg)

# ruleorder: lastdb_index > lastal_alignment

rule all: # run all rules ######################################################
    input:
        expand(
            "{dataDir}/{{sample}}/{{sample}}.msa".format(dataDir=dataDir),
            sample = samples
        ),
        expand(
            "{dataDir}/{{sample}}/checkSpecifity/allOligos_reverse.fasta".format(dataDir=dataDir),
            sample = samples
        )
        # expand(
        #     "{dataDir}/{{sample}}/{{sample}}.index.{{ext}}".format(dataDir=dataDir),
        #     sample = samples,
        #     ext = ["bck", "des", "prj", "sds", "ssp", "suf", "tis"]
        # ),
        # expand(
        #     "{dataDir}/{{sample}}/checkSpecifity/seal.out".format(dataDir=dataDir),
        #     sample = samples
        # ),
        # expand(
        #     "{dataDir}/{{sample}}/checkSpecifity/seal.out.stats".format(dataDir=dataDir),
        #     sample = samples
        # ),

include: "rules/all_preprocessing.smk"
include: "rules/reverse_kmer_dsk.smk"
include: "rules/reverse_oligo_filter.smk"
include: "rules/reverse_check_specifity.smk"
