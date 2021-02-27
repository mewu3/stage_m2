#!/usr/bin/env python3.8
import pandas as pd

configfile: "config.yaml"
# samples = pd.read_table(config["samples"], header=0)
# samples_dir = config["samples_dir"]
# results_dir = config["results_dir"]
samples=["enterovirus", "coronavirus"]

include: "rules/remove-duplicate-seq.smk"
include: "rules/multiple-seq-alignment.smk"
include: "rules/split_overlap_chunks.smk"
# include: "rules/kmer-counting.smk"
