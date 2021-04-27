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
        f"{dataDir}/{{sample}}/kmer{kmerSize}/evaluation/allOligos_reverse.set.coverage"
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
        dict_idPosi = defaultdict(list)
        # dict_idSpecieCount = defaultdict(lambda: defaultdict(int))
        dict_aceIdSpecie = defaultdict(str)
        dict_speciesCount = defaultdict(int)

        with open(input[0], "r") as file:
            for li in file.readlines()[1:]:
                li = li.rstrip("\n").split()
                oligoId.append(f"p{li[0]}")
                position = f"{li[-2]}-{li[-1]}"
                dict_idPosi[li[0]] = li[8:]

        df_bowtie = pd.read_table(input[2], comment="@", sep="\t", names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "OPT1", "OPT2", "OPT3", "OPT4"])

        df_bowtie = df_bowtie[df_bowtie["QNAME"].isin(oligoId)][["QNAME", "RNAME"]]

        dict_idAce = {k: g["RNAME"].tolist() for k, g in df_bowtie.groupby("QNAME")}

        for id in dict_idAce:
            dict_idAce[id] = list(map(lambda x: x.split(".")[0], dict_idAce[id]))

        input2_open = open(input[1], "r")
        for li in input2_open.readlines()[1:]:
            li = li.rstrip("\n").split("\t")
            aceId = li[1]
            specie = li[4]
            dict_aceIdSpecie[aceId]=specie
            dict_speciesCount[specie] += 1
        input2_open.close()

        outputOpen = open(output[0], "w")
        for id in dict_idAce:
            dict = defaultdict(int)
            for ace in  dict_idAce[id]:
                specie = dict_aceIdSpecie[ace]
                dict[specie]+=1
            id = id.lstrip("p")
            posi = dict_idPosi[id]
            posi = "\t".join(map(str, dict_idPosi[id]))
            for specie in dict:
                species_count = dict[specie]
                totolCount = dict_speciesCount[specie]
                outputOpen.write(f"{posi}\t{specie}\t{species_count}\t{totolCount}\n")
        outputOpen.close()
