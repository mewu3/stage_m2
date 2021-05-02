rule evaluation1:
    input:
        f"{dataDir}/{{sample}}/{{sample}}.uniq" if config["curated"] != "true" else lambda wildcards: config["samples"][wildcards.sample],
        file_aceIDtaxID,
        file_taxIDLineage
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
        df_taxLineage = dd.read_csv(input[2], sep=",", header=0, assume_missing=True, usecols=["tax_id", "species", "no rank"], dtype={"tax_id":"int64", "species":"object", "no rank":"object"})
        df_taxID = df_taxID[df_taxID["accession"].isin(acessionID)]
        df_out = df_taxID.merge(df_taxLineage, how="left", left_on="taxid", right_on="tax_id")

        df_out.compute().to_csv(output[0], sep="\t")

rule evaluation2:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allOligo.set",
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/aceID-taxID-species.tsv"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/evaluation/allOligo.set.coverage"
    params:
        splitFilesDir = f"{dataDir}/{{sample}}/splitFiles"
    run:
        import os
        import sys
        import re
        from collections import defaultdict
        from Bio import SeqIO
        from Bio.Seq import Seq

        input1 = input[0]
        input2 = input[1]
        output1 = output[0]

        dict_posiSeq = defaultdict(str)
        dict_posiAce = defaultdict(list)
        dict_aceIdSpecie = defaultdict(str)
        dict_speciesCount = defaultdict(int)

        input1_open = open(input1, "r")
        for line in input1_open.readlines()[1:]:
            line = line.rstrip("\n").split()
            position = f"{line[-2]}-{line[-1]}"
            oligo = line[1]
            dict_posiSeq[position] = oligo
        input1_open.close()

        for posi in dict_posiSeq:
            seq = str(Seq(dict_posiSeq[posi]).reverse_complement())
            splitFiles = f"{params.splitFilesDir}/reverse{posi}.fasta"
            accessionID_ls = []
            record_dict = SeqIO.to_dict(SeqIO.parse(splitFiles, "fasta"))
            for record in record_dict:
                if re.search(seq, str(record_dict[record].seq), re.I):
                    acessionID = record_dict[record].id.split(".")[0]
                    accessionID_ls.append(acessionID)
            dict_posiAce[posi]=accessionID_ls

        input2_open = open(input2, "r")
        for li in input2_open.readlines()[1:]:
            li = li.rstrip("\n").split("\t")
            aceId = li[1]
            specie = li[4]
            dict_aceIdSpecie[aceId]=specie
            dict_speciesCount[specie] += 1
        input2_open.close()

        output1_open = open(output1, "w")
        for posi in dict_posiAce:
            dict=defaultdict(int)
            for specie in dict_speciesCount:
                dict[specie]=0
            for acessionID in dict_posiAce[posi]:
                spec = dict_aceIdSpecie[acessionID]
                dict[spec]+=1
            for specie in dict:
                species_count = dict[specie]
                totolCount = dict_speciesCount[specie]
                posi = "\t".join(map(str, posi.split("-")))
                output1_open.write(f"{posi}\t{specie}\t{species_count}\t{totolCount}\n")
        output1_open.close()

rule evaluation3:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allKmercount.sorted.calculated.txt",
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allKmerCount.sorted.calculated.filtered.txt",
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allKmerCount.sorted.calculated.filtered.spec.txt",
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allOligo.set"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/evaluation/allOligo_before.tsv",
        f"{dataDir}/{{sample}}/kmer{kmerSize}/evaluation/allOligo_after1.tsv",
        f"{dataDir}/{{sample}}/kmer{kmerSize}/evaluation/allOligo_after2.tsv",
        f"{dataDir}/{{sample}}/kmer{kmerSize}/evaluation/allOligo_after3.tsv"
    run:
        import sys
        import pandas as pd
        from collections import defaultdict

        for index in range(len(input)):
            df = pd.read_table(input[index], sep="\t", header=0)
            df = df.sort_values(["start", "kmerCount"], ascending=False)
            df_first = df.groupby("start", as_index=False).first()
            df_first = df_first[["oligo", 'kmerCount', 'CG%', 'Entropy', 'Tm', 'homodimer-dG', 'hairpin-dG', 'start', 'end']]
            df_first.to_csv(output[index] , sep="\t", index=False)
