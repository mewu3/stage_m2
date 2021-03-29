#!/usr/bin/env python3.8
import os

configfile: "config.yaml"

samples = list(config["samples"].keys())
dataDir = config["dataDir"]
refSeq = config["refSeq"]

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


# def aggregate_input(wildcards):
#     checkpoint_output = checkpoints.splitFiles.get(**wildcards).output[0]
#     return expand(f"{dataDir}/{{sample}}/filtering/reverse{{seg}}.calculated2.filtered",
#                   sample = wildcards.sample,
#                   seg = glob_wildcards(os.path.join(checkpoint_output, "reverse{seg}.fasta")).seg)


rule all:
    input:
        # expand(
        #     f"{dataDir}/{{sample}}/checkSpecifity/references.DB.{{ext}}",
        #     sample = samples,
        #     ext = ["nhr", "nin", "nog", "nsd", "nsi", "nsq"]
        #
        # )
        expand(
            f"{dataDir}/{{sample}}/checkSpecifity/allOligos_reverse.filtered.spec.fasta",
            sample = samples
        )

include: "rules/all_preprocessing.smk"
include: "rules/reverse_kmer_dsk.smk"
include: "rules/reverse_kmer_filter.smk"
include: "rules/reverse_kmer_spec.smk"
