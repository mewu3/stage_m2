#!/usr/bin/env python3.8
from snakemake.utils import validate, min_version
import pandas as pd

### define minimum required version for snakmake
# min_version("5.32.2")

### load config file and access sample file through "sample.tsv"
configfile: "config.yaml"
# validate(config, schema="") # validate whether config file is adequate
# sample = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(sample, schema="")

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

rule
