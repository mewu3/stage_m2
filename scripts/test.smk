#!/usr/bin/env python3.8
import os

configfile: "../config.yaml"
samples = list(config["samples"].keys())
file_refSeq = config["file_refSeq"]
file_aceIDtaxID = config["file_aceIDtaxID"]
file_taxIDLineage = config["file_taxIDLineage"]
dataDir = config["dataDir"]
kmerSize = config["kmerSize"]


### msspe basic 
# rule all: 
#     input: 
#         expand(
#             f"{dataDir}/{{sample}}/kmer{kmerSize}/evaluation/allKmerCount.sorted.calculated.filtered.allOligo.set.coverage",
#             sample = samples 
#         )

# rule test:
#     input:
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/evaluation/allOligo_after2.tsv",
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/aceID-taxID-species.tsv"
#     output:
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/evaluation/allKmerCount.sorted.calculated.filtered.allOligo.set.coverage"
#     params:
#         splitFilesDir = f"{dataDir}/{{sample}}/splitFiles"
#     run:
#         import os
#         import sys
#         import re
#         from collections import defaultdict
#         from Bio import SeqIO
#         from Bio.Seq import Seq

#         input1 = input[0]
#         input2 = input[1]
#         output1 = output[0]

#         dict_posiSeq = defaultdict(str)
#         dict_posiAce = defaultdict(list)
#         dict_aceIdSpecie = defaultdict(str)
#         dict_speciesCount = defaultdict(int)

#         input1_open = open(input1, "r")
#         for line in input1_open.readlines()[1:]:
#             line = line.rstrip("\n").split()
#             position = f"{line[-2]}-{line[-1]}"
#             oligo = line[0]
#             dict_posiSeq[position] = oligo
#         input1_open.close()

#         for posi in dict_posiSeq:
#             seq = str(Seq(dict_posiSeq[posi]).reverse_complement())
#             splitFiles = f"{params.splitFilesDir}/reverse{posi}.fasta"
#             accessionID_ls = []
#             record_dict = SeqIO.to_dict(SeqIO.parse(splitFiles, "fasta"))
#             for record in record_dict:
#                 if re.search(seq, str(record_dict[record].seq), re.I):
#                     acessionID = record_dict[record].id.split(".")[0]
#                     accessionID_ls.append(acessionID)
#             dict_posiAce[posi]=accessionID_ls

#         input2_open = open(input2, "r")
#         for li in input2_open.readlines()[1:]:
#             li = li.rstrip("\n").split("\t")
#             aceId = li[1]
#             specie = li[4]
#             dict_aceIdSpecie[aceId]=specie
#             dict_speciesCount[specie] += 1
#         input2_open.close()

#         output1_open = open(output1, "w")
#         for posi in dict_posiAce:
#             dict=defaultdict(int)
#             for specie in dict_speciesCount:
#                 dict[specie]=0
#             for acessionID in dict_posiAce[posi]:
#                 spec = dict_aceIdSpecie[acessionID]
#                 dict[spec]+=1
#             for specie in dict:
#                 species_count = dict[specie]
#                 totolCount = dict_speciesCount[specie]
#                 posi = "\t".join(map(str, posi.split("-")))
#                 output1_open.write(f"{posi}\t{specie}\t{species_count}\t{totolCount}\n")
#         output1_open.close()

### msspe variant 
rule all: 


rule getKmerPosition1:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allKmerCount.sorted.calculated.filtered.txt"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.fasta"
    run:
        import os
        from Bio.Seq import Seq

        outFile = open(output[0], "w")
        input1Open = open(input[0], "r")
        for li in input1Open.readlines()[1:]:
            li = li.rstrip("\n")
            ls = li.split()
            id = ls[0]
            oligo = str(Seq(ls[1]).reverse_complement())
            outFile.write(f">p{id}\n{oligo}\n")
        input1Open.close()
        outFile.close()

rule getKmerPosition2:
    input:
        f"{dataDir}/{{sample}}/{{sample}}{clusterIdentity}.msa"
    output:
        multiext(
            f"{dataDir}/{{sample}}/{{sample}}{clusterIdentity}.msa.DB.", "1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"
        )
    params:
        out = f"{dataDir}/{{sample}}/{{sample}}{clusterIdentity}.msa.DB",
        threads = config["thread"]
    shell:
        """
        bowtie-build \
        --threads {params.threads} \
        {input} {params.out}
        """

rule getKmerPosition3:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.fasta",
        lambda wildcards: expand(
            f"{dataDir}/{{sample}}/{{sample}}{clusterIdentity}.msa.DB.{{ext}}",
            ext = ["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"],
            sample = wildcards.sample
        )
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.bowtie"
    params:
        refDB = f"{dataDir}/{{sample}}/{{sample}}{clusterIdentity}.msa.DB",
        threads = config["thread"]
    shell:
        """
        bowtie \
        -x {params.refDB} -p {params.threads} \
        -v 0 -l 7 --norc -a --sam \
        -f {input[0]} {output}
        """

rule getKmerPosition4:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.spec.bowtie",
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allKmerCount.sorted.calculated.filtered.spec.txt"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.spec.position"
    run:
        import os
        import pandas as pd
        from collections import defaultdict

        dict_idPosiCount = defaultdict(lambda: defaultdict(int))
        dict_idPosi = defaultdict(int)

        df_bowtie = pd.read_table(input[0], comment="@", sep="\t", names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "OPT1", "OPT2", "OPT3", "OPT4"])

        Match = df_bowtie[df_bowtie["FLAG"] == 0]["QNAME"].tolist()
        Match = set(map(lambda x : int(x.lstrip("p")), Match))

        position = df_bowtie[df_bowtie["FLAG"] == 0].groupby(["QNAME", "POS"]).size().reset_index(name = "POScount")
        position = position.sort_values(["QNAME","POScount"], ascending=False).groupby("QNAME").first()
        position.index = position.index.str.replace("p", "").astype("int64")

        df = pd.read_table(input[1], sep="\t", header=0, index_col=0)
        df_filtered = df[df.index.isin(Match)]
        df_out = df_filtered.join(position)
        df_out.to_csv(output[0], sep='\t', index=True)

rule checkHeterodimer:
    input:
        f"{dataDir}/{{sample}}/{{sample}}{clusterIdentity}.msa",
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.spec.position"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allOligo.set"
    params:
        step = config["step"],
        kmerSize = config["kmerSize"],
        deltaG = config["heterodimer-deltaG"],
        monovalentConc = config["dimer-oligotm"]["monovalent-conc"],
        divalentConc = config["dimer-oligotm"]["divalent-conc"],
        dNTPConc = config["dimer-oligotm"]["dNTP-conc"],
        dnaConc = config["dimer-oligotm"]["dna-conc"],
        thermodynamicPara = config["dimer-oligotm"]["thermodynamic-para"],
        saltCorrelation = config["dimer-oligotm"]["salt-correlation"]
    run:
        import pandas as pd
        from Bio import SeqIO
        import numpy as np
        import primer3
        from collections import defaultdict

        step = params.step
        kmerSize = params.kmerSize

        def calculate_heterodimer(kmer1, kmer2):
            heterodimer = primer3.calcHeterodimer(kmer1, kmer2,
                                                  mv_conc = params.monovalentConc,
                                                  dv_conc = params.divalentConc,
                                                  dntp_conc = params.dNTPConc,
                                                  dna_conc = params.dnaConc).dg
            if float(heterodimer) > params.deltaG: # superior than -9000 - no 2nd structure formation
                return True


        records = list(SeqIO.parse(input[0], "fasta"))
        seqLength = len(records[0].seq)

        chunks = [[i,i+step] for i in range(0, seqLength, step)]

        for chunk in chunks:
            if chunk[-1] > seqLength:
                chunk[-1] = seqLength

        df = pd.read_table(input[1], sep="\t", header=0, index_col=0)
        df.index.name = "index"
        df = df.drop("POScount", 1)
        df = df.rename(columns={"POS":"start"})
        df["end"] = df["start"] + kmerSize -1
        df["chunk-start"] = ""
        df["chunk-end"] = ""

        criteria =[]
        for chunk in chunks:
            criteria.append(df.end.between(chunk[0], chunk[1]))

        chunk_start = [x[0] for x in chunks]
        chunk_end = [x[1] for x in chunks]

        df["chunk-start"] = np.select(criteria, chunk_start, 0)
        df["chunk-end"] = np.select(criteria, chunk_end, 0)

        df = df.sort_values(["kmerCount"], ascending=False).groupby("chunk-end").head(3)

        dict_idInfo = df.to_dict("index")
        df["index"] = df.index
        # dict_posiIDs = {k: g["index"].tolist() for k, g in df[["index", "chunk-end"]].groupby("chunk-end")} # lose kmerCount order

        dict_posiIDs=defaultdict(list)

        for index in df.index:
            position = df["chunk-end"][index]
            dict_posiIDs[position].append(index)

        heterodimer=[]
        dict_heterodimer = defaultdict(list)

        for p1 in dict_posiIDs:
            for id1 in dict_posiIDs[p1]:
                kmer1 = dict_idInfo[id1]["oligo"]
                for p2 in dict_posiIDs:
                    if p1 != p2:
                        for id2 in dict_posiIDs[p2]:
                            kmer2 = dict_idInfo[id2]["oligo"]
                            if calculate_heterodimer(kmer1, kmer2):
                                heterodimer.append([p1, id1, id2])
                                dict_heterodimer[id1].append(id2)

        oligoSet = []
        seenPosition = []
        oligoSet.append(heterodimer[0][1])
        seenPosition.append(heterodimer[0][0])

        for hetero in heterodimer:
            posi = hetero[0]
            id = hetero[1]
            if posi not in seenPosition:
                seenPosition.append(posi)
                if id not in oligoSet:
                    for register in oligoSet:
                        if id in dict_heterodimer[register]:
                            if id not in oligoSet:
                                oligoSet.append(id)

        df_filtered = df[df.index.isin(oligoSet)].drop("index", axis=1)
        df_filtered = df_filtered.sort_values("chunk-end")
        df_filtered.to_csv(output[0], sep='\t', index=True)
