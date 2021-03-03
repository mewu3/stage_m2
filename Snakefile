#!/usr/bin/env python3.8
import glob

configfile: "config.yaml"

samples = list(config["samples"].keys())
datadir = config["datadir"]

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

def aggregate_input(wildcards):
    checkpoint_out = checkpoints.split_overlap_chunks.get(**wildcards).output[0]


rule all: # run all rules ######################################################
    input:
        ### remove duplicate sequences from input fasta
        expand(
            datadir+"/duplicate_removed/{sample}.uniq.fasta",
            sample = samples
        ),
        ### generate MSA files in format fasta
        expand(
            datadir + "/msa/{sample}.msa.fasta",
               sample = samples
        ),
        ### generate split files
        dynamic(
            expand(
                datadir + "/splitFiles/{sample}.forward.{seg}.fasta",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/splitFiles/{sample}.reverse.{seg}.fasta",
                sample = samples,
                seg="{seg}"
            )
        ),
        aggregate_input

include: "rules/remove_duplicate.smk"
include: "rules/multiple_seq_alignment.smk"
include: "rules/splitFiles.smk"
include: "rules/kmer_dsk.smk"
# include: "rules/kmer_kmc.smk"
# include: "rules/oligo_filter.smk"
