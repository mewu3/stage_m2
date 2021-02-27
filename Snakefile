#!/usr/bin/env python3.8
import pandas as pd

configfile: "config.yaml"
# samples = pd.read_table(config["samples"], header=0)
# samples_dir = config["samples_dir"]
# results_dir = config["results_dir"]
samples = ["enterovirus", "coronavirus"]
datapath = config["data"]

def get_fasta():
    return "{}/*.fasta".format(datapath)

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


rule seqkit: # remove duplicate sequences ######################################
    input:
        # "test/{sample}.fasta"
        get_fasta()
    output:
        "{datapath}/results/duplicate_removed/{sample}.uniq.fasta"
    conda:
        "envs/seqkit.yaml"
    log:
        "{datapath}/results/duplicate_removed/seqkit.{sample}.log"
    shell:
        "cat {input} |seqkit rmdup -s -o {output}"



rule mafft: # multiple sequences alignment #####################################
    input:
        "test/results/duplicate_removed/{sample}.uniq.fasta"
    output:
        "test/results/msa/{sample}.msa.fasta"
    log:
        "test/results/msa/mafft.{sample}.log"
    conda:
        "envs/mafft.yaml"
    params:
        threads = config["mafft"]["threads"],
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
        "mafft \
        --thread {params.threads} \
        {params.algorithm} \
        --op {params.op} \
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

rule split_overlap_chunks: # split MSA fasta file into seperates files #########
    input:
        "test/results/msa/{sample}.msa.fasta"
    output:
        dir = directory("test/results/split_files/")
    conda:
        "env/pyfaidx.yaml"
    params:
        step = config["step"]
        overlap = config["overlap"]
    run:
        import os
        from pyfaidx import Fasta

        fasta_file = input[0]
        step = params.step
        overlap = params.overlap

        os.makedirs(output.dir)

        seqLength = len(fasta_file[0])

        remainder = seqLength % (step-overlap)
        chunkNumber = int(seqLength / (step-overlap))
        print("The sequences are splitted into {} chunks, and there are {} bp left".format())

rule dsk:
    input:
        dir = "test/splitFiles/"
    output:
        dir = directory("test/kmerCounting/dsk")
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
        histo = config["dsk"]["histo"]
    script:
        "scripts/dsk.py"

rule dskOutput:
    input:
        dir = "test/kmerCounting/dsk"
    output:
        dir = directory("test/kmerCounting/dsk_output")
    script:
        "scripts/dsk_output.py"

rule kmc:
    input:
        dir = "test/splitFiles"
    output:
        dir = directory("test/kmerCounting/kmc3")
    # log:
    #     "test/kmerCounting/kmc3/log.txt"
    conda:
        "envs/kmc3.yaml"
    params:
        kmerSize = config["kmc3"]["kmer-size"]
    script:
        "scripts/kmc3.py"

rule kmcOutput:
    input:
        dir = "test/kmerCounting/kmc3"
    output:
        dir = directory("test/kmerCounting/kmc3_output")
    script:
        "scripts/kmc3_output.py"
    script:
        "scripts/dsk.py"
