rule clustering:
    input:
        f"{dataDir}/{{sample}}/{{sample}}.uniq"
    output:
        f"{dataDir}/{{sample}}/{{sample}}{clusterIdentity}.uniq"
    params:
        identity = config["cd-hit"]["identity"],
        threads = config["thread"],
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
        f"{dataDir}/{{sample}}/{{sample}}{clusterIdentity}.uniq"
    output:
        f"{dataDir}/{{sample}}/{{sample}}{clusterIdentity}.msa"
    log:
        f"{dataDir}/{{sample}}/log/mafft.log"
    conda:
        "envs/mafft.yaml"
    params:
        threads = config["thread"],
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

rule getKmerPosition1:
    input:
        # f"{dataDir}/{{sample}}/{{kmerSize}}/allKmerCount.sorted.calculated.filtered.spec.txt"
        f"{dataDir}/{{sample}}/{{kmerSize}}/allKmerCount.sorted.calculated.txt"
    output:
        # f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.spec.fasta"
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.fasta"
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
        # f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.spec.fasta",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.fasta",
        lambda wildcards: expand(
            f"{dataDir}/{{sample}}/{{sample}}{clusterIdentity}.msa.DB.{{ext}}",
            ext = ["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"],
            sample = wildcards.sample
        )
    output:
        # f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.spec.bowtie"
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.bowtie"
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
        # f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.spec.bowtie",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.bowtie",
        # f"{dataDir}/{{sample}}/{{kmerSize}}/allKmerCount.sorted.calculated.filtered.spec.txt"
        f"{dataDir}/{{sample}}/{{kmerSize}}/allKmerCount.sorted.calculated.txt"
    output:
        #f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.spec.position"
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.position"
    run:
        import os
        import pandas as pd
        import numpy as np
        from collections import defaultdict

        dict_idPosiCount = defaultdict(lambda: defaultdict(int))
        dict_idPosi = defaultdict(int)

        df_bowtie = pd.read_table(input[0], comment="@", sep="\t", names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "OPT1", "OPT2", "OPT3", "OPT4"])

        # Match = df_bowtie[df_bowtie["FLAG"] == 0]["QNAME"].tolist()
        # Match = set(map(lambda x : int(x.lstrip("p")), Match))

        position = df_bowtie[df_bowtie["FLAG"] == 0].groupby(["QNAME", "POS"]).size().reset_index(name = "POScount")
        position = position.sort_values(["QNAME","POScount"], ascending=False).groupby("QNAME").first()
        position.index = position.index.str.replace("p", "").astype("int64")
        position=position[["POS", "POScount"]]

        df = pd.read_table(input[1], sep="\t", header=0, index_col=0)
        # df_filtered = df[df.index.isin(Match)]
        df_out = df.merge(position, how="left", left_index=True, right_index=True)
        # df_out.fillna(0)
        # df_out["POS"] = df_out["POS"].replace(np.nan, 0)
        # print(df_out)
        df_out.to_csv(output[0], sep='\t', index=True)

rule checkHeterodimer:
    input:
        f"{dataDir}/{{sample}}/{{sample}}{clusterIdentity}.msa",
        # f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.spec.position"
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.position",
        f"{dataDir}/{{sample}}/{{kmerSize}}/allKmerCount.sorted.calculated.filtered.spec.txt"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/allOligo.set"
    params:
        step = config["step"],
        kmerSize = lambda wildcards: config["kmerSize"][wildcards.kmerSize],
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

        df_before = pd.read_table(input[1], sep="\t", header=0, index_col=0)
        df_before.index.name = "index"
        df_before = df_before.drop("POScount", 1)
        df_before = df_before.rename(columns={"POS":"start"})
        df_before["end"] = df_before["start"] + kmerSize -1
        df_before["chunk-start"] = ""
        df_before["chunk-end"] = ""

        criteria =[]
        for chunk in chunks:
            criteria.append(df_before.end.between(chunk[0], chunk[1]))

        chunk_start = [x[0] for x in chunks]
        chunk_end = [x[1] for x in chunks]

        df_before["chunk-start"] = np.select(criteria, chunk_start, 0)
        df_before["chunk-end"] = np.select(criteria, chunk_end, 0)

        df_position = df_before[["start", "end", "chunk-start", "chunk-end"]]

        df_after2 = pd.read_csv(input[2], sep="\t", header=0, index_col=0)
        df_after2 = df_after2.merge(df_position, how="left", left_index=True, right_index=True)

        df = df_after2.sort_values(["kmerCount"], ascending=False).groupby("chunk-end").head(3)

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
        df_filtered.fillna(0, inplace=True)
        df_filtered.to_csv(output[0], sep='\t', index=True)
