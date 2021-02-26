#!/usr/bin/env python3.8
import pandas as pd

# load config file and access sample file through "sample.tsv" #################
configfile: "config.yaml"
sample = pd.read_table(config["samples"]).set_index("target", drop=False)

# SNAKEMAKE RULE ###############################################################

## MMSEQ2: remove duplicate sequences
rule mmseq2: # cluster input sequences
    input:
        "test/enterovirus_50.fasta"
    output:
        "test/enterovirus_50_clsutered.fasta"
    log:
        "test/mmseq2.log"
    shell:
        "mmseqs createdb {input} | \
         mmseq2 easy-linclust {output} "

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
        "test/mafft/sequences_aligned.fasta"
    log:
        "test/mafft/mafft.log"
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

### Split aligment fasta file on overlapping chunks
rule split_overlap_chunks:
    input:
        fasta = "test/mafft/sequences_aligned.fasta"
    output:
        dir = directory("test/splitFiles/")
    script:
        "scripts/overlap_segs.py"

### Kmer counting
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
        histoMax = config["dsk"]["histo-max"],# with os.scandir("test/splitFiles") as it:
#     for entry in it:
#         if not entry.name.startswith('.') and entry.is_file():
#             print(entry.name)
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
