#!/usr/bin/env python3.8
import os

configfile: "config.yaml"

samples = list(config["samples"].keys())
kmerSize = list(config["kmerSize"].keys())

file_refSeq = config["file_refSeq"]
file_aceIDtaxID = config["file_aceIDtaxID"]
file_taxIDLineage = config["file_taxIDLineage"]

dataDir = config["dataDir"]

deduplication = config["deduplication"]
clustering = config["clustering"]

clusterIdentity = config["cd-hit"]["identity"]

rule all:
    input:
        "/datater/wu/data/genomeVariability/enterovirus_fastaANI"


rule evaluation1:
    input:
        "/datater/wu/data/enterovirusCurated.fasta",
        file_aceIDtaxID,
        file_taxIDLineage
    output:
        "/datater/wu/data/genomeVariability/"
    run:
        import os
        import pandas as pd
        import dask.dataframe as dd
        import numpy as np
        import time

        acessionID = os.popen(f"cat {input[0]}|grep '^>'|cut -f1 -d' '|sed 's/>//g'").read().split("\n")
        acessionID = list(map(lambda x : x.split(".")[0], acessionID))

        df_taxID = dd.read_csv(input[1], sep="\t", header=0, usecols=["accession", "taxid"], dtype={"accession":"object", "taxid":"int64"})
        df_taxLineage = dd.read_csv(input[2], sep=",", header=0, assume_missing=True, usecols=["tax_id", "species", "no rank"], dtype={"tax_id":"int64", "species":"object", "no rank":"object"})
        df_taxID = df_taxID[df_taxID["accession"].isin(acessionID)]
        df_out = df_taxID.merge(df_taxLineage, how="left", left_on="taxid", right_on="tax_id")

        df_out.compute().to_csv(output[0], sep="\t")

rule test:
    input:
        "/datater/wu/data/genomeVariability/enterovirus.txt"
    output:
        "/datater/wu/data/genomeVariability/enterovirus_fastaANI"
    shell:
        "lib/FastANI/fastANI --rl {input} --ql {input} -o {output} --matrix -t 12"
