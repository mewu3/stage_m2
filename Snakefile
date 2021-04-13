#!/usr/bin/env python3.8
import os

configfile: "config.yaml"

samples = list(config["samples"].keys())
dataDir = config["dataDir"]
refSeq = config["refSeq"]
kmerSize = config["jellyfish-count"]["kmer-size"]

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
        expand(
            f"{dataDir}/{{sample}}/{{sample}}.uniq",
            sample = samples
        ),
        expand(
            f"{dataDir}/{{sample}}/{{sample}}.msa",
            sample = samples
        ),
        expand(
            f"{dataDir}/{{sample}}/filtering{kmerSize}/allOligos_reverse.filtered",
            sample = samples
        ),
        expand(
            f"{dataDir}/{{sample}}/checkSpecifity{kmerSize}/allOligos_reverse.filtered.fasta",
            sample = samples
        ),
        expand(
            f"{dataDir}/{{sample}}/checkSpecifity{kmerSize}/references.DB.{{ext}}",
            sample = samples,
            ext = ["nhr", "nin", "nog", "nsd", "nsi", "nsq"]

        ),
        expand(
            f"{dataDir}/{{sample}}/checkSpecifity{kmerSize}/allOligos_reverse.filtered.blastn.out",
            sample = samples
        ),
        expand(
            f"{dataDir}/{{sample}}/checkSpecifity{kmerSize}/allOligos_reverse.filtered.spec.fasta",
            sample = samples
        ),
        expand(
            f"{dataDir}/{{sample}}/dimer{kmerSize}/allOligos_reverse.filtered.spec.heterodimer",
            sample = samples
        ),
        expand(
            f"{dataDir}/{{sample}}/dimer{kmerSize}/allOligos_reverse.set",
            sample = samples
        ),
        expand(
            f"{dataDir}/{{sample}}/dimer{kmerSize}/allOligos_reverse.set.fasta",
            sample = samples
        ),
        expand(
            f"{dataDir}/{{sample}}/evaluation{kmerSize}/aceID-taxID-species.tsv",
            sample = samples
        ),
        expand(
            f"{dataDir}/{{sample}}/evaluation{kmerSize}/allOligos_reverse.set.coverage",
            sample = samples
        ),
        expand(
            f"{dataDir}/{{sample}}/evaluation{kmerSize}/allOligo_before.tsv",
            sample = samples
        ),
        expand(
            f"{dataDir}/{{sample}}/evaluation{kmerSize}/allOligo_after1.tsv",
            sample = samples
        ),
        expand(
            f"{dataDir}/{{sample}}/evaluation{kmerSize}/allOligo_after2.tsv",
            sample = samples
        ),
        expand(
            f"{dataDir}/{{sample}}/evaluation{kmerSize}/allOligo_after3.tsv",
            sample = samples
        )


include: "rules/all_preprocessing.smk"
include: "rules/reverse_kmer_kmerCount.smk"
include: "rules/reverse_kmer_filter.smk"
include: "rules/reverse_kmer_spec.smk"
include: "rules/reverse_kmer_dimer.smk"
include: "rules/evaluation.smk"
