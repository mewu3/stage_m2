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
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.clustered.spec.position"
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
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.clustered.spec.position"
    output: 
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/test"
    params: 
        step = config["step"]
    run: 
        import os
        import pandas as pd
        import numpy as np 
        from Bio import SeqIO
        from collections import defaultdict

        dict_id_length = {}

        df = pd.read_table(input[2], comment="@", sep="\t", names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "OPT1", "OPT2", "OPT3", "OPT4"])
        df = df[["QNAME", "RNAME", "POS"]]
        
        dict_records = SeqIO.index(input[0], "fasta")
        ls_length = []
        for id in dict_records: 
            dict_id_length[id] = len(dict_records[id])
            ls_length.append(len(dict_records[id]))
        
        df["REFLENGTH"] = df["RNAME"].map(dict_id_length)
        df["POSITION"] = df["POS"]/df["REFLENGTH"] * 100 
        df = df[["QNAME", "RNAME", "POSITION"]]

        meanLen = int(sum(ls_length) / len(ls_length))

        chunksNumber = len([[i,i+params.step] for i in range(0, meanLen, params.step)])      
        step = round(100/chunksNumber, 1)
        chunks = [[i,i+step] for i in np.arange(0, 100, step)]
        print(chunks)

        # for chunk in chunks:
        #     if chunk[-1] > meanLen:
        #         chunk[-1] = meanLen




        # dict_id_position = {k: g["POS"].tolist() for k, g in df.groupby("QNAME")} 
        # print(dict_id_position)

