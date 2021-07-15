rule kmerSelecting1:
    input:
        input_kmerCounting1
    output:
        multiext(
            f"{dataDir}/{{sample}}/{{sample}}.DB.", "1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"
        )
    params:
        out = f"{dataDir}/{{sample}}/{{sample}}.DB",
        threads = config["thread"]
    shell:
        """
        bowtie-build \
        --threads {params.threads} \
        {input[0]} {params.out}
        """

rule kmerSelecting2:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.clustered.spec.fasta",
        lambda wildcards: expand(
            f"{dataDir}/{{sample}}/{{sample}}.DB.{{ext}}",
            sample = wildcards.sample,
            ext = ["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"]
        )
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.clustered.spec.bowtie"
    params:
        refDB = f"{dataDir}/{{sample}}/{{sample}}.DB",
        threads = config["thread"]
    shell:
        """
        bowtie \
        -x {params.refDB} -f {input[0]} \
        -v 0 -l 7 -a --nofw \
        --sam \
        -p {params.threads} \
        {output}
        """

rule kmerSelecting3:
    input:
        input_kmerCounting1,
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.clustered.spec.fasta",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.clustered.spec.bowtie",
        f"{dataDir}/{{sample}}/{{kmerSize}}/allKmerCount.sorted.calculated.filtered.txt"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/allOligo.set",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allOligo.set.position"
    params:
        step = config["step"],
        deltaG = config["heterodimer-deltaG"],
        monovalentConc = config["dimer-oligotm"]["monovalent-conc"],
        divalentConc = config["dimer-oligotm"]["divalent-conc"],
        dNTPConc = config["dimer-oligotm"]["dNTP-conc"],
        dnaConc = config["dimer-oligotm"]["dna-conc"],
        thermodynamicPara = config["dimer-oligotm"]["thermodynamic-para"],
        saltCorrelation = config["dimer-oligotm"]["salt-correlation"]
    run:
        import os
        import math
        import pandas as pd
        import numpy as np
        from Bio import SeqIO
        import primer3
        from collections import defaultdict
        import time

        step = params.step

        ls_heterodimer = []
        ls_id_toExclude = []
        ls_oligo = set()

        dict_id_count = {}
        dict_aceID_length = {}
        dict_aceID_seq = {}
        dict_posi_ids = {}

        def calculate_heterodimer(kmer1, kmer2):
            heterodimer = primer3.calcHeterodimer(kmer1, kmer2,
                                                  mv_conc = params.monovalentConc,
                                                  dv_conc = params.divalentConc,
                                                  dntp_conc = params.dNTPConc,
                                                  dna_conc = params.dnaConc).dg
            if float(heterodimer) > params.deltaG: # superior than -9000 - no 2nd structure formation
                return True

        ls_heterodimer = []
        dict_id_oligo = SeqIO.index(input[1], "fasta")
        for id1 in dict_id_oligo:
            seq1 = str(dict_id_oligo[id1].seq)
            for id2 in dict_id_oligo:
                seq2 = str(dict_id_oligo[id2].seq)
                if not calculate_heterodimer(seq1, seq2):
                    heterodimer = [id1, id2]
                    ls_heterodimer.append(heterodimer)

        ### remove heterodimer incompaticility ### 
        with open(input[3], "r") as fi:
            for l in fi.readlines()[1:]:
                l=l.rstrip("\n").split()
                dict_id_count[l[0]] = l[2]

        for heterodimer in ls_heterodimer:
            id1 = heterodimer[0]
            id2 = heterodimer[1]
            count1 = dict_id_count[id1]
            count2 = dict_id_count[id2]
            if count1 > count2:
                ls_id_toExclude.append(count2)
            else:
                ls_id_toExclude.append(count1)

        dict_aceID_seq = SeqIO.index(input[0], "fasta")
        for aceID in dict_aceID_seq:
            dict_aceID_length[aceID] = len(dict_aceID_seq[aceID])

        df = pd.read_table(input[2], comment="@", sep="\t", names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "OPT1", "OPT2", "OPT3", "OPT4"])
        df = df[["QNAME", "RNAME", "POS"]]
        df["REFLENGTH"] = df["RNAME"].map(dict_aceID_length)
        df["POSITION"] = df["POS"]/df["REFLENGTH"] * 100
        df = df[["QNAME", "RNAME", "POSITION"]]
        df = df[~df.index.isin(ls_id_toExclude)]

        ### inference roughly position for oligonucleotides
        meanLen = int(sum(list(dict_aceID_length.values())) / len(list(dict_aceID_length.values())))

        chunks = [[i,i+step] for i in range(0, meanLen, step)]
        chunks[-1][-1] = meanLen
        chunkNumber = len(chunks)
        remainder = chunks[-1][-1]-chunks[-1][0]

        if remainder <= step/2 :
            step = int(remainder/chunkNumber) +1 + step

        chunks = [[i,i+step] for i in range(0, meanLen, step)]
        chunks[-1][-1] = meanLen
        chunkNumber = len(chunks)

        step = math.floor(100/chunkNumber)
        chunks = [[i,i+step] for i in np.arange(0, 100, step)]
        chunks[-1][-1] = 100
        remainder = chunks[-1][-1]-chunks[-1][0]

        if remainder <= step/2 :
            step = int(remainder/chunkNumber) +1 + step

        chunks = [[i,i+step] for i in range(0, 100, step)]
        chunks[-1][-1] = 100
        chunkNumber = len(chunks)

        criteria = []
        for chunk in chunks:
            criteria.append(df.POSITION.between(chunk[0], chunk[1]))

        chunk_start = [x[0] for x in chunks]
        chunk_end = [x[1] for x in chunks]

        df["START"] = ""
        df["START"]=np.select(criteria, chunk_start, 0)
        df["END"] = ""
        df["END"]=np.select(criteria, chunk_end, 0)

        df["COUNT"] = df.groupby("QNAME")["QNAME"].transform('count')
        df = df[["QNAME", "RNAME", "START", "END", "COUNT"]]

        ### select oligonucleotides ### 
        for start, group_df in df.groupby("START"):

            group_df = group_df.sort_values(["COUNT", "QNAME"], ascending=False)

            nb_refs = len(group_df["RNAME"].unique())

            count = 0
            while count < nb_refs*0.9:
                index_max = group_df["COUNT"].idxmax()
                id_max = group_df.loc[index_max].QNAME
                ls_oligo.add(id_max)
                ls_refs = group_df[group_df["QNAME"] == id_max]["RNAME"].to_list()
                group_df = group_df[~group_df.RNAME.isin(ls_refs)]
                count += len(ls_refs)

        ls_oligo = list(ls_oligo)

        df_ = df[df.QNAME.isin(ls_oligo)]
        df_.to_csv(output[1], sep="\t", index=False)

        df_oligo = pd.read_table(input[3], sep="\t", header=0, index_col=0)
        df_oligo = df_oligo[df_oligo.index.isin(ls_oligo)]
        df_oligo.to_csv(output[0], sep="\t", index=True)
