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

rule all: # run all rules ######################################################
    input:
        # remove duplicate sequences
        expand(
            datadir+"/duplicate_removed/{sample}.uniq.fasta",
            sample = samples
        ),
        # multiple sequences alignment
        expand(
            datadir + "/msa/{sample}.msa.fasta",
               sample = samples
        ),
        # split msa fasta
        dynamic(
            expand(
                datadir + "/splitFiles/{sample}/forward/{seg}.fasta",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/splitFiles/{sample}/reverse/{seg}.fasta",
                sample = samples,
                seg="{seg}"
            )
        ),
        ### forward
        dynamic(
            expand(
                datadir + "/kmerCounting/{sample}/forward/{seg}.h5",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/kmerCounting/{sample}/forward/{seg}_dsk.txt",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/kmerCounting/{sample}/forward/{seg}_dsk.sort",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/filtering/{sample}/forward/{seg}_dsk.table",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/filtering/{sample}/forward/{seg}_dsk.tsv",
                sample = samples,
                seg="{seg}"
            )
        ),
        ### reverse
        dynamic(
            expand(
                datadir + "/kmerCounting/{sample}/reverse/{seg}.h5",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/kmerCounting/{sample}/reverse/{seg}_dsk.txt",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/kmerCounting/{sample}/reverse/{seg}_dsk.sort",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/filtering/{sample}/reverse/{seg}_dsk.table",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/filtering/{sample}/reverse/{seg}_dsk.tsv",
                sample = samples,
                seg="{seg}"
            )
        )

include: "rules/remove_duplicate.smk"
include: "rules/multiple_seq_alignment.smk"
include: "rules/splitFiles.smk"
include: "rules/kmer_dsk.smk"
# include: "rules/kmer_kmc.smk"
include: "rules/oligo_filter.smk"
