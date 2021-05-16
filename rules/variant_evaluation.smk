rule evaluation1:
    input:
        # lambda wildcards: config["samples"][wildcards.sample] if config["deduplication"] else f"{dataDir}/{{sample}}/{{sample}}.uniq",
        f"{dataDir}/{{sample}}/{{sample}}{clusterIdentity}.msa",
        file_aceIDtaxID,
        file_taxIDLineage
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/aceID-taxID-species.tsv"
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
        f"{dataDir}/{{sample}}/{{kmerSize}}/allOligo.set",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/aceID-taxID-species.tsv",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.bowtie"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/evaluation/allOligo.set.coverage"
    run:
        import os
        import sys
        import re
        from collections import defaultdict
        from Bio import SeqIO
        from Bio.Seq import Seq
        import pandas as pd
        import numpy as np

        oligoId = []
        dict_idPosi = defaultdict(list)
        # dict_idSpecieCount = defaultdict(lambda: defaultdict(int))
        dict_aceIdSpecie = defaultdict(str)
        dict_speciesCount = defaultdict(int)

        with open(input[0], "r") as file:
            for li in file.readlines()[1:]:
                li = li.rstrip("\n").split()
                oligoId.append(f"p{li[0]}")
                # position = f"{li[-2]}-{li[-1]}"
                dict_idPosi[li[0]] = li[10:]

        df_bowtie = pd.read_table(input[2], comment="@", sep="\t", names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "OPT1", "OPT2", "OPT3", "OPT4"])
        df_bowtie = df_bowtie[df_bowtie["FLAG"] != 4 ]
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
            for specie in dict_speciesCount:
                dict[specie]=0
            for ace in dict_idAce[id]:
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

rule evaluation3:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/allOligo.set"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/evaluation/allOligo_after3.tsv"
    shell:
        "cp {input[0]} {output[0]}"

rule evaluation4:
    input:
        f"{dataDir}/{{sample}}/{{sample}}{clusterIdentity}.msa",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.position",
        f"{dataDir}/{{sample}}/{{kmerSize}}/allKmerCount.sorted.calculated.filtered.txt",
        f"{dataDir}/{{sample}}/{{kmerSize}}/allKmerCount.sorted.calculated.filtered.spec.txt"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/evaluation/allOligo_before.tsv",
        f"{dataDir}/{{sample}}/{{kmerSize}}/evaluation/allOligo_after1.tsv",
        f"{dataDir}/{{sample}}/{{kmerSize}}/evaluation/allOligo_after2.tsv"
    params:
        step = config["step"],
        kmerSize = lambda wildcards: config["kmerSize"][wildcards.kmerSize],
    run:
        import pandas as pd
        from Bio import SeqIO
        import numpy as np
        from collections import defaultdict

        step = params.step
        kmerSize = params.kmerSize

        records = list(SeqIO.parse(input[0], "fasta"))
        seqLength = len(records[0].seq)

        chunks = [[i,i+step] for i in range(0, seqLength, step)]

        for chunk in chunks:
            if chunk[-1] > seqLength:
                chunk[-1] = seqLength

        df_before = pd.read_table(input[1], sep="\t", header=0, index_col=0)
        df_before.index.name = "index"
        df_before = df_before.drop("POScount", 1)
        df_before = df_before.rename(columns={"POS":"rstart"})
        df_before["rend"] = df_before["rstart"] + kmerSize -1
        df_before["start"] = ""
        df_before["end"] = ""

        criteria =[]
        for chunk in chunks:
            criteria.append(df_before.rend.between(chunk[0], chunk[1]))

        chunk_start = [x[0] for x in chunks]
        chunk_end = [x[1] for x in chunks]

        df_before["start"] = np.select(criteria, chunk_start, 0)
        df_before["end"] = np.select(criteria, chunk_end, 0)

        df_position = df_before[["rstart", "rend", "start", "end"]]

        df_before = df_before.sort_values(["kmerCount"], ascending=False).groupby("end", as_index=False).head(1)


        df_after1 = pd.read_csv(input[2], sep="\t", header=0, index_col=0)
        df_after1 = df_after1.merge(df_position, how="left", left_index=True, right_index=True)
        df_after1 = df_after1.sort_values(["kmerCount"], ascending=False).groupby("end", as_index=False).head(1)

        df_after2 = pd.read_csv(input[3], sep="\t", header=0, index_col=0)
        df_after2 = df_after2.merge(df_position, how="left", left_index=True, right_index=True)
        df_after2 = df_after2.sort_values(["kmerCount"], ascending=False).groupby("end", as_index=False).head(1)

        df_before.fillna(0, inplace=True)
        df_after1.fillna(0, inplace=True)
        df_after2.fillna(0, inplace=True)

        df_before.to_csv(output[0], sep="\t", index=True)
        df_after1.to_csv(output[1], sep="\t", index=True)
        df_after2.to_csv(output[2], sep="\t", index=True)
