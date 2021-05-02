rule checkHeterodimer:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allKmerCount.sorted.calculated.filtered.spec.txt"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allOligo.set"
    params:
        monovalentConc = config["dimer-oligotm"]["monovalent-conc"],
        divalentConc = config["dimer-oligotm"]["divalent-conc"],
        dNTPConc = config["dimer-oligotm"]["dNTP-conc"],
        dnaConc = config["dimer-oligotm"]["dna-conc"],
        thermodynamicPara = config["dimer-oligotm"]["thermodynamic-para"],
        saltCorrelation = config["dimer-oligotm"]["salt-correlation"],
        deltaG = config["heterodimer-deltaG"]
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
            if float(heterodimer) > params.deltaG: # superior than -9000 - no 2nd structure formation
                return True

        df = pd.read_table(input[0], sep="\t", header=0, index_col=0)

        df = df.sort_values(["kmerCount"], ascending=False).groupby("end").head(3)

        dict_idInfo = df.to_dict("index")

        df["index"] = df.index

        dict_posiIDs=defaultdict(list)

        for index in df.index:
            position = df["end"][index]
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
        df_filtered = df_filtered.sort_values("end")
        df_filtered.to_csv(output[0], sep='\t', index=True)
