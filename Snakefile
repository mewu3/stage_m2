#!/usr/bin/env python3.8
from snakemake.utils import validate, min_version
import pandas as pd

# define minimum required version for snakmake #################################
min_version("5.32.2")

# load config file and access sample file through "sample.tsv" #################
configfile: "config.yaml"
# validate(config, schema="") # validate whether config file is adequate
# sample = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(sample, schema="")

# environment variables ########################################################
envvars:

# SNAKEMAKE RULE ###############################################################
### MMSEQ2
# rule mmseq2: # cluster input sequences
#     input:
#         "test/enterovirus_50.fasta"
#     output:
#         "test/enterovirus_50_clsutered.fasta"
#     log:
#         "test/mmseq2.log"
#     shell:
#         "mmseqs createdb {input} | \
#          mmseq2 easy-linclust {output} "

### MAFFT: multiple sequences alignment
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

rule mafft:
    input:
        "test/sequences.fasta"
    output:
        "test/sequences_aligned.fasta"
    log:
        "test/mafft.log"
    conda:
        "envs/mafft.yaml"
    params:
        algorithm = config["mafft"]["algorithm"],
        op = config["mafft"]["op"],
        ep = config["mafft"]["ep"],
        bl = config["mafft"]["bl"],
        jtt = config["mafft"]["jtt"],
        tm = config["mafft"]["tm"],
        fmodel = config["mafft"]["fmodel"],
        clustalout = config["mafft"]["clustalout"],
        inputorder = config["mafft"]["inputorder"],
        reorder = config["mafft"]["reorder"],
        treeout = config["mafft"]["treeout"],
        quiet = config["mafft"]["quiet"]
    shell:
        "mafft {params.algorithm} \
        --op {params.opvous} \
        --ep {params.ep} \
        --bl {params.bl} \
        --jtt {params.jtt} \
        --tm {params.tm} \
        {params.fmodel} \
        {params.clustalout} \
        {params.inputorder} \
        {params.reorder} \
        {params.treeout} \
        {params.quiet} \
        {input} > {output} \
        2> {log}"

rule split_overlap_chunks:
    input:
        fasta = "test/sequences_aligned.fasta"
    output:
        dir = directory("test/splitFiles/")
    script:
        "scripts/overlap_segs.py"

rule count_kmers:
    input:
        dir = directory("test/splitFiles/")
    output:
        dir = directory("test/kmerCounting/dsk")
    log:
        "test/dsk_output"
    params:
        nbCores = config["dsk"]["nb-cores"],
        maxMemory = config["dsk"]["max-memory"],
        maxDisk = config["dsk"]["max-disk"],
        outCompress = config["dsk"]["out-compress"],
        storage = config["dsk"]["storage-type"],
        verbose = config["dsk"]["verbose"],
        kmerSize = config["dsk"]["kmer-size"],
        abundanceMin = config["dsk"]["abundance-min"],
        abundanceMax = config["dsk"]["abundance-max"],
        abundanceMinThreshold = config["dsk"]["abundance-min-threshold"],
        solidityKind = config["dsk"]["solidity-kind"],
        solidityCustom = config["dsk"]["solidity-custom"],
        solideKmerOut = config["dsk"]["solid-kmers-out"],
        histoMax = config["dsk"]["histo-max"],
        histo2D = config["dsk"]["histo2D"],
        histo = config["dsk"]["histo"],

    # run:
    #     import os
    #     import sys
    #
    #     outputDir = os.makedirs(output.dir)
    #
    #     with os.scandir(input[0]) as it:
    #         for entry in it:
    #             if not entry.name.startswith('.') and entry.is_file():
    #                 os.system("""lib/dsk/build/bin/dsk
    #                           -nb-cores {params.nbCores}
    #                           -max-memory {params.maxMemory}
    #                           -max-disk {params.maxDisk}
    #                           -out-compress {params.outCompress}
    #                           -storage-type {params.storage}
    #                           -verbose {params.verbose}
    #                           -kmer-size {params.kmerSize}
    #                           -abundance-min {params.abundanceMin}
    #                           -abundance-max {params.abundanceMax} 
    #                           -abundance-min-threshold {params.abundanceMinThreshold}
    #                           -solidity-kind {params.solidityKind}
    #                           -solidity-custom {params.solidityCustom}
    #                           -solid-kmers-out {params.solideKmerOut}
    #                           -histo-max {params.histoMax}
    #                           -histo2D {params.histo2D}
    #                           -histo {params.histo}
    #                           -file""" + entry.path + "-out-dir" + outputDir )

    shell: 
