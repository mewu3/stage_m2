from snakemake.utils import validate, min_version
import pandas as pd
# from conf import load_sample_table

### define minimum required version for snakmake
min_version("5.32.2")

### load config file and access sample file data through "samples.tsv"
# configfile: "config.yaml"
# validate(config, schema="") # validate whether config file is adequate
# samples = load_sample_table
# validate(sample, schema="") 

### snakemake rules
rule multiple-sequences-alignment:
    input:
    output:
    conda:
        "envs/mafft.yaml"
    shell:
