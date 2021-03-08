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
        ### remove duplicate sequences
        expand(
            datadir+"/duplicate_removed/{sample}.uniq.fasta",
            sample = samples
        ),
        ### multiple sequences alignment
        expand(
            datadir + "/msa/{sample}.msa.fasta",
               sample = samples
        ),
        ### forward
        dynamic(
            expand(
                datadir + "/{sample}/split_forward/{seg}.fasta",
                sample = samples,
                seg="{seg}"
            )
        ),
        ### dsk kmer counting
        dynamic(
            expand(
                datadir + "/{sample}/dsk/forward_{seg}.h5",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/{sample}/dsk/forward_{seg}.kCount",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/{sample}/dsk/forward_{seg}.kCountSorted",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/{sample}/primer3_Tm/forward_{seg}.Tm",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/{sample}/primer3_Tm/forward_{seg}.TmFilterered",
                sample = samples,
                seg="{seg}"
            )
        ),
        expand(
            datadir + "/{sample}/primer3_Tm/all_oligo",
            sample = samples
        ),
        expand(
            datadir + "/{sample}/primer3_Tm/all_oligo.fasta",
            sample = samples
        ),
        expand(
            datadir + "/{sample}/primer3_Tm/all_oligo.fasta.fai"
        )

include: "rules/remove_duplicate.smk"
include: "rules/multiple_seq_alignment.smk"
include: "rules/splitFiles.smk"

include: "rules/kmer_dsk.smk"
include: "rules/dsk_oligo_filter.smk"

# include: "rules/kmer_kmc.smk"
