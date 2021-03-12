#!/usr/bin/env python3.8
import glob

configfile: "config.yaml"

samples = list(config["samples"].keys())
datadir = config["datadir"]
kmerCounter = "dsk"

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
        ### remove duplicate sequences
        expand(
            datadir+"/duplicate_removed/{sample}.uniq",
            sample = samples
        ),
        ### multiple sequences alignment
        expand(
            datadir + "/msa/{sample}.msa",
               sample = samples
        ),
        ### forward
        # dynamic(
        #     expand(
        #         datadir + "/{sample}/splitFiles/forward{seg}.fasta",
        #         sample = samples,
        #         seg="{seg}"
        #     )
        # ),
        dynamic(
            expand(
                datadir + "/{sample}/splitFiles/reverse{seg}.fasta",
                sample = samples,
                seg="{seg}"
            )
        ),
        ### dsk kmer counting reverse
        dynamic(
            expand(
                datadir + "/{sample}/dsk/reverse{seg}.h5",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/{sample}/dsk/reverse{seg}.kCount",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/{sample}/dsk/reverse{seg}.kCountSorted",
                sample = samples,
                seg="{seg}"
            )
        ),
        ### dsk kmer counting forward
        # dynamic(
        #     expand(
        #         datadir + "/{sample}/dsk/forward{seg}.h5",
        #         sample = samples,
        #         seg="{seg}"
        #     )
        # ),
        # dynamic(
        #     expand(
        #         datadir + "/{sample}/dsk/forward{seg}.kCount",
        #         sample = samples,
        #         seg="{seg}"
        #     )
        # ),
        # dynamic(
        #     expand(
        #         datadir + "/{sample}/dsk/forward{seg}.kCountSorted",
        #         sample = samples,
        #         seg="{seg}"
        #     )
        # ),
        ### reverse calculating parameters
        dynamic(
            expand(
                datadir + "/{sample}/filtering/reverse{seg}.calculated",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/{sample}/filtering/reverse{seg}.fa",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/{sample}/filtering/reverse{seg}.nessieOut",
                sample = samples,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                datadir + "/{sample}/filtering/reverse{seg}.calculated2",
                sample = samples,
                seg="{seg}"
            )
        ),
        expand(
            datadir + "/{sample}/filtering/allOligos_reverse.calculated2",
            sample = samples,
        ),
        expand(
            datadir + "/{sample}/filtering/allOligos_reverse.filtered",
            sample = samples,
            seg="{seg}"
        ),
        expand(
            datadir + "/{sample}/checkSpecifity/allOligos_reverse.fasta",
            sample = samples
        ),
        # expand(
        #     datadir + "/{sample}/checkSpecifity/all_oligo.fasta",
        #     sample = samples
        # ),
        # expand(datadir + "/{sample}/primer3_Tm/all_oligo.fasta.{ext}",
        #     sample = samples,
        #     ext = ["fai", "json", "log", "primerqc", "primerqc.fai"]
        # ),
        # expand(
        #     datadir + "/{sample}/hairpin_out",
        #     sample = samples
        # ),
        # expand(
        #     datadir + "/{sample}/primer3_ntthal/allOligo_hairpin",
        #     sample = samples
        # ),
        # expand(
        #     datadir + "{sample}/ptrimer3_ntthal/allOligo_hairpin_dimer",
        #     sample = samples
        # )

include: "rules/all_remove_duplicate.smk"
include: "rules/all_multiple_seq_alignment.smk"
include: "rules/all_splitFiles.smk"
include: "rules/reverse_kmer_dsk.smk"
include: "rules/reverse_oligo_filter.smk"
include: "rules/reverse_check_specifity.smk"
