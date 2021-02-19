from snakemake.utils import validate, min_version
import pandas as pd
# from conf import load_sample_table

### define minimum required version for snakmake
min_version("5.32.2")

### load config file and access sample file through "sample.tsv"
configfile: "config.yaml"
# validate(config, schema="") # validate whether config file is adequate
# sample = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(sample, schema="")



# SNAKEMAKE RULE ###############################################################
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


if config["mafft"]["clustalout"] == "on":
    config["mafft"]["clustalout"] = "--clustalout"
if config["mafft"]["inputorder"] == "on":
    config["mafft"]["inputorder"] = "--inputorder"
if config["mafft"]["reorder"] == "on":
    config["mafft"]["reorder"] = "--reorder"
if config["mafft"]["treeout"] == "on":
    config["mafft"]["treeout"] = "--treeout"
if config["mafft"]["quiet"] == "on":
    config["mafft"]["quiet"] = "--quiet"

rule mafft: # multiple sequence alignment
    input:
        "test/sequences.fasta"
    output:
        "test/sequences_aligned.clustal"
    log:
        "test/mafft.log"
    conda:
        "envs/mafft.yaml"
    params:
        algorithm = config["mafft"]["algorithm"],
        outformat = config["mafft"]["clustalout"],
        scorematrix = config ["mafft"]["tm"]
    shell:
        "mafft {params.algorithm} \
        --tm {params.scorematrix} \
        {params.outformat} \
        {input} > {output} \
        2> {log}"
