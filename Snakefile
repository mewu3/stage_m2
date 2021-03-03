#!/usr/bin/env python3.8
import glob

configfile: "config.yaml"

samples = list(config["samples"].keys())
datadir = config["datadir"]
kmerCounter = "kmc"

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

kmerCounting_prefix = [os.path.splitext(f)[0] for f in os.listdir(datadir + "/splitFiles") if f.endswith(".fasta")]

rule all: # run all rules ######################################################
    input:
        # remove duplicate sequences from input fasta
        expand(
            datadir+"/duplicate_removed/{sample}.uniq.fasta",
            sample = samples
        ),
        # generate MSA files in format fasta
        expand(
            datadir + "/msa/{sample}.msa.fasta",
               sample = samples
        ),
        # generate split files
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
        # kmerCounting with dsk
        expand(
            datadir + "/kmerCounting/dsk/{prefix}.h5",
            prefix = kmerCounting_prefix
        ),
        expand(
            datadir + "/kmerCounting/dsk/{prefix}.txt",
            prefix = kmerCounting_prefix
        ),
        expand(
            datadir + "/kmerCounting/dsk/{prefix}.sort",
            prefix = kmerCounting_prefix
        ),
        # kmerCounting with kmc3
        expand(
            datadir + "/kmerCounting/kmc/{prefix}.kmc_suf",
            prefix = kmerCounting_prefix
        ),
        expand(
            datadir + "/kmerCounting/kmc/{prefix}.kmc_pre",
            prefix = kmerCounting_prefix
        ),
        expand(
            datadir + "/kmerCounting/kmc/{prefix}.txt",
            prefix = kmerCounting_prefix
        ),
        expand(
            datadir + "/kmerCounting/kmc/{prefix}.sort",
            prefix = kmerCounting_prefix
        ),
        expand(
            datadir + "/filter/kmc/{prefix}.txt",
            prefix = kmerCounting_prefix
        ),
        expand(
            datadir + "/filter/" + kmerCounter + "/{prefix}.table",
            prefix = kmerCounting_prefix
        )

include: "rules/remove_duplicate.smk"
include: "rules/multiple_seq_alignment.smk"
include: "rules/splitFiles.smk"
include: "rules/kmer_dsk.smk"
include: "rules/kmer_kmc.smk"
include: "rules/oligo_filter.smk"
