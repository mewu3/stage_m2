rule evaluation1:
    input:
        f"{dataDir}/{{sample}}/{{sample}}.uniq",
        taxID,
        taxLineage
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/aceID-taxID-species.tsv"
    run:
        import os
        import pandas as pd
        import dask.dataframe as dd
        import numpy as np
        import time

        acessionID = os.popen(f"cat {input[0]}|grep '^>'|cut -f1 -d' '|sed 's/>//g'").read().split("\n")
        acessionID = list(map(lambda x : x.split(".")[0], acessionID))

        df_taxID = dd.read_csv(input[1], sep="\t", header=0, usecols=["accession", "taxid"], dtype={"accession":"object", "taxid":"int64"})
        df_taxLineage = dd.read_csv(input[2], sep=",", header=0, assume_missing=True, usecols=["tax_id", "species"], dtype={"tax_id":"int64", "species":"object"})
        df_taxID = df_taxID[df_taxID["accession"].isin(acessionID)]
        df_out = df_taxID.merge(df_taxLineage, how="left", left_on="taxid", right_on="tax_id")

        df_out.compute().to_csv(output[0], sep="\t")

rule evaluation2:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allOligo.set",
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/aceID-taxID-species.tsv",
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.spec.bowtie"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allOligos_reverse.set.coverage"
    params:
        splitFilesDir = f"{dataDir}/{{sample}}/splitFiles"
    run:
        import os
        import sys
        import re
        from collections import defaultdict
        from Bio import SeqIO
        from Bio.Seq import Seq
        import pandas as pd

        oligoId = []

        with open(input[0], "r") as file:
            for li in file.readlines()[1:]:
                li = li.rstrip("\n").split()
                oligoId.append(f"p{li[0]}")

        df_bowtie = pd.read_table(input[2], comment="@", sep="\t", names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "OPT1", "OPT2", "OPT3", "OPT4"])

        print(df_bowtie[df_bowtie["QNAME"].isin(oligoId)][["QNAME", "RNAME"]].to_dict)
