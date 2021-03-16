#!/usr/bin/env python3.8
import os

configfile: "config.yaml"

dataDir = config["dataDir"]
sampleDir = config["sampleDir"]

sample = config["sample"]
kmerCounter = config["kmerCounter"]

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
#     return expand(dataDir + "/" + sample + "/splitFiles/forward{seg}.calculated2",
#                   seg = glob_wildcards(os.path.join(checkpoint_output, "forward{seg}.calculated2")).seg)
#
# def aggregate_reverseInput(wildcards):
#     checkpoint_output = checkpoints.splitFiles.get(**wildcards).output[0]
#     return expand(dataDir + "/" + sample + "/splitFiles/reverse{seg}.calculated2",
#                   seg = glob_wildcards(os.path.join(checkpoint_output, "reverse{seg}.calculated2")).seg)

# def aggregate_reverseInput(wildcards):
#     checkpoint_output = checkpoints.splitFiles.get(**wildcards).output[0]
#     return expand("{dataDir}/{sample}/splitFiles/reverse{{seg}}.calculated2".format(dataDir = dataDir, sample = sample),
#                   seg = glob_wildcards(os.path.join(checkpoint_output, "reverse{seg}.calculated2")).seg)

rule all: # run all rules ######################################################
    input:
        # dynamic(
        #     expand(
        #         "{dataDir}/{sample}/filtering/reverse{{seg}}.calculated2".format(dataDir = dataDir,
        #                                                                          sample = sample),
        #         seg = "{seg}"
        #     )
        # ),
        "{}/{}/filtering/allOligos_reverse.calculated2".format(dataDir, sample)

include: "rules/all_preprocessing.smk"
include: "rules/reverse_kmer_dsk.smk"
include: "rules/reverse_oligo_filter.smk"
# include: "rules/reverse_check_specifity.smk"
