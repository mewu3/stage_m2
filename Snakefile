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

rule all:
    input: 
        expand(
            f"{dataDir}/{{sample}}/filtering/allOligos_reverse.filtered", 
            sample = samples
        )

include: "rules/all_preprocessing.smk"
include: "rules/reverse_kmer_dsk.smk"
include: "rules/reverse_oligo_filter.smk"
include: "rules/reverse_check_specifity.smk"
