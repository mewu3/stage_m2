rule clustering:
    input:
        f"{dataDir}/{{sample}}/{{sample}}.uniq"
    output:
        f"{dataDir}/{{sample}}/{{sample}}{clusterIden}.uniq"
    params:
        identity = config["cd-hit"]["identity"],
        threads = config["cd-hit"]["thread"],
        memory = config["cd-hit"]["memory"]
    shell:
        "./lib/cd-hit-v4.8.1-2019-0228/cd-hit-est \
        -i {input} \
        -o {output} \
        -c {params.identity} \
        -T {params.threads} \
        -M {params.memory}"

rule MSA:
    input:
        f"{dataDir}/{{sample}}/{{sample}}{clusterIden}.uniq"
    output:
        f"{dataDir}/{{sample}}/{{sample}}{clusterIden}.msa"
    log:
        f"{dataDir}/{{sample}}/log/mafft.log"
    conda:
        "envs/mafft.yaml"
    params:
        threads = config["mafft"]["threads"],
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
        --thread {params.threads} \
        --op {params.op} \
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

# rule getKmerPosition:
#     input:
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/allKmerCount.sorted.calculated.filtered.spec.txt",
#         f"{dataDir}/{{sample}}/{{sample}}{clusterIden}.msa"
#     output:
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.spec.position"
#     params:
#         step = config["step"],
#         overlap = config["overlap"]
#     run:
#         from Bio import SeqIO
#         from itertools import islice
#         import re
#         from collections import defaultdict
#         from Bio.Seq import Seq
#         import pandas as pd
#
#         step = params.step
#         overlap = params.overlap
#
#         dict_idKmer = defaultdict(str)
#         dict_posiIdCount = defaultdict(lambda: defaultdict(int))
#
#         input1Open = open(input[0], "r")
#         for li in input1Open.readlines()[1:]:
#             li = li.rstrip("\n")
#             ls = li.split()
#             id = ls[0]
#             kmer = ls[1]
#             dict_idKmer[id]=kmer
#         input1Open.close()
#
#         record_dict = SeqIO.to_dict(SeqIO.parse(input[1], "fasta"))
#         for record in record_dict:
#             for id in dict_idKmer:
#                 kmer = Seq(dict_idKmer[id])
#                 kmer = str(kmer.reverse_complement())
#                 iterator = re.finditer(kmer, str(record_dict[record].seq), re.I)
#                 for match in iterator:
#                     position = match.span()
#                     dict_posiIdCount[position][id]+=1
#                     print(id, position)
#
#         key1 = list(record_dict)[0]
#         length = len(record_dict[key1].seq)
#
#         chunks = [[i,i+step] for i in range(0, length, step)]
#
#         df = pd.DataFrame(columns=["start", "end", "ID", "count", "chunk-start", "chunk-end"])
#         df_original = pd.read_table(input[0], sep="\t", header=0,index_col=0)
#
#         count = 0
#         for posi in dict_posiIdCount:
#             count += 1
#             start, end = posi
#             start = int(start)
#             end = int(end)
#             for id in dict_posiIdCount[posi]:
#                 c = dict_posiIdCount[posi][id]
#                 df.loc[count] = pd.Series({"start":start, "end":end, "ID":id, "count":c, "chunk-sart":"unknow", "chunk-end": "unknow"})
#
#         for chunk in chunks:
#             if chunk[-1] > length:
#                 chunk[-1] = length
#
#         for index in df.index:
#             start = df["start"][index]
#             end = df["end"][index]
#             for chunk in chunks:
#                 rg = range(chunk[0], chunk[1])
#                 if end in rg:
#                     df["chunk-start"][index] = chunk[0]
#                     df["chunk-end"][index] = chunk[1]
#
#         df["ID"] = df["ID"].astype("int64")
#         df = df.set_index("ID")
#
#         df_output = df.join(df_original)
#         # df_output = df_output.sort_values(["chunk-end", "count"], ascending=False)
#         df_output = df_output.sort_values(["count", "chunk-end"], ascending=False)
#         df_out = df_output.groupby("chunk-end").head(3)
#
#         df_out.to_csv(output[0], sep="\t", index=True)
#
# rule checkHeterodimer:
#     input:
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.spec.position",
#     output:
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/allOligo.set"
#     params:
#         monovalentConc = config["oligotm"]["monovalent-conc"],
#         divalentConc = config["oligotm"]["divalent-conc"],
#         dNTPConc = config["oligotm"]["dNTP-conc"],
#         dnaConc = config["oligotm"]["dna-conc"],
#         thermodynamicPara = config["oligotm"]["thermodynamic-para"],
#         saltCorrelation = config["oligotm"]["salt-correlation"]
#     run:
#         from collections import defaultdict
#         import primer3
#
#         dict_IDInfo = defaultdict(list)
#         dict_posiIDs = defaultdict(list)
#         dict_heterodimer = defaultdict(list)
#
#         heterodimer = []
#         oligoSet = []
#         seenPosition=[]
#
#         def calculate_heterodimer(kmer1, kmer2):
#             heterodimer = primer3.calcHeterodimer(kmer1, kmer2,
#                                                   mv_conc = params.monovalentConc,
#                                                   dv_conc = params.divalentConc,
#                                                   dntp_conc = params.dNTPConc,
#                                                   dna_conc = params.dnaConc).dg
#             if float(heterodimer) > -9000: # superior than -9000 - no 2nd structure formation
#                 return True
#
#         input1Open = open(input[0], "r")
#         for li in input1Open.readlines()[1:]:
#             li = li.rstrip("\n")
#             ls = li.split()
#             id = ls[0]
#             start = ls[-2]
#             end = int(start) +
#             chunkStart = ls[4]
#             chunkEnd = ls[5]
#             kmer = ls[6]
#             kmerCount = ls[7]
#             dict_IDInfo[id]=[kmer, kmerCount, start, end, chunkStart, chunkEnd]
#             dict_posiIDs[chunkEnd].append(id)
#         input1Open.close()
#
#         output1Open = open(output[0], "w")
#         output1Open.write(f"Parent Child\n")
#
#         for p1 in dict_posiIDs:
#             # seenPosition.append(p1)
#             for id1 in dict_posiIDs[p1]:
#                 kmer1 = dict_IDInfo[id1][0]
#                 for p2 in dict_posiIDs:
#                     if p1 != p2:
#                         for id2 in dict_posiIDs[p2]:
#                             kmer2 = dict_IDInfo[id2][0]
#                             if calculate_heterodimer(kmer1, kmer2):
#                                 heterodimer.append([id1, id2])
#                                 dict_heterodimer[id1].append(id2)
#                                 output1Open.write(f"{id1} {id2}\n")
#
#         output1Open.close()
#
#         input1Open = open(input[0], "r")
#         oligoSet.append(input1Open.readlines()[1].split()[0])
#         input1Open.close()
#
#         input1Open = open(input[0], "r")
#         for li in input1Open.readlines()[1:]:
#             li = li.rstrip("\n")
#             ls = li.split()
#             id = ls[0]
#             posi = ls[5]
#             if posi not in seenPosition:
#                 seenPosition.append(posi)
#                 if id not in oligoSet:
#                     for register in oligoSet:
#                         if id in dict_heterodimer[register]:
#                             if id not in oligoSet:
#                                 oligoSet.append(id)
#                                 # print(id, register)
#                                 # print(oligoSet)
#         input1Open.close()
#
#         output1Open=open(output[0], "w")
#         for id in oligoSet:
#             ls = dict_IDInfo[id]
#             formatedList = "\t".join(map(str,ls))
#             output1Open.write(f"{id}\t{formatedList}\n")
#         output1Open.close()

rule getKmerPosition1:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allKmerCount.sorted.calculated.filtered.spec.txt"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.spec.fasta"
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
        f"{dataDir}/{{sample}}/{{sample}}{clusterIden}.msa"
    output:
        multiext(
            f"{dataDir}/{{sample}}/{{sample}}{clusterIden}.msa.DB.", "1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"
        )
    params:
        out = f"{dataDir}/{{sample}}/{{sample}}{clusterIden}.msa.DB",
        threads = 12
    shell:
        """
        bowtie-build \
        --threads {params.threads} \
        {input} {params.out}
        """

rule getKmerPosition3:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.spec.fasta"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.spec.bowtie"
    params:
        refDB = f"{dataDir}/{{sample}}/{{sample}}{clusterIden}.msa.DB",
        threads = 12
    shell:
        """
        bowtie \
        -x {params.refDB} -f {input} \
        -v 0 -l 7 --ff \
        -a \
        --sam \
        -p {params.threads} \
        {output}
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

        # print(pd.read_table(input[0], comment="@"))
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
        f"{dataDir}/{{sample}}/{{sample}}{clusterIden}.msa",
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.spec.position"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allOligo.set"
    params:
        step = config["step"],
        kmerSize = config["jellyfish"]["kmer-size"],
        monovalentConc = config["oligotm"]["monovalent-conc"],
        divalentConc = config["oligotm"]["divalent-conc"],
        dNTPConc = config["oligotm"]["dNTP-conc"],
        dnaConc = config["oligotm"]["dna-conc"],
        thermodynamicPara = config["oligotm"]["thermodynamic-para"],
        saltCorrelation = config["oligotm"]["salt-correlation"]
    run:
        import pandas as pd
        from Bio import SeqIO
        import numpy as np
        import primer3
        from collections import defaultdict

        def calculate_heterodimer(kmer1, kmer2):
            heterodimer = primer3.calcHeterodimer(kmer1, kmer2,
                                                  mv_conc = params.monovalentConc,
                                                  dv_conc = params.divalentConc,
                                                  dntp_conc = params.dNTPConc,
                                                  dna_conc = params.dnaConc).dg
            if float(heterodimer) > -9000: # superior than -9000 - no 2nd structure formation
                return True

        step = params.step
        kmerSize = params.kmerSize

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

        # outputOpen = open(output[0], "w")
        # outputOpen.write(f"ID\toligo\tkmerCount\tCG%\tEntropy\tTm\thomodier-dG\thairpin-dG\tstart\tend\tchunk-start\tchunk-end\n")
        # for id in oligoSet:
        #     ls = [v for k, v in dict_idInfo[id].items()]
        #     formatedList = "\t".join(map(str,ls))
        #     outputOpen.write(f"{id}\t{formatedList}\n")
        #
        # outputOpen.write(f"# The aligned genome length is {seqLength}.\n")
        # outputOpen.close()
