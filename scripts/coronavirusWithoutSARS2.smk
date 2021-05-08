# !/usr/bin/env python3.8
import os

sars2Ref = "/datater/wu/data/tmp/sars2_variants.fasta"
datadir = ["MSSPE-variant", "MSSPE-basic"]
sample = ["coronavirusWithoutSarsCov2", "coronavirusCuratedWithoutSarsCov2"]
kmerSize = ["kmer15", "kmer13"]

configfile: "config.yaml"

def fetchConfigParameters():
    if config["mafft"]["fmodel"] == "on":
        config["mafft"]["fmodel"] = "--fmodel"
    else:
        config["mafft"]["fmodel"] = ""

    if config["mafft"]["clustalout"] == "on":
        config["mafft"]["clustalout"] = "--clustalout"
    else:
        config["mafft"]["clustalout"] = ""

    if config["mafft"]["inputorder"] == "on":
        config["mafft"]["inputorder"] = "--inputorder"
    else :
        config["mafft"]["inputorder"] = ""

    if config["mafft"]["reorder"] == "on":
        config["mafft"]["reorder"] = "--reorder"
    else:
        config["mafft"]["reorder"] = ""

    if config["mafft"]["treeout"] == "on":
        config["mafft"]["treeout"] = "--treeout"
    else :
        config["mafft"]["treeout"] = ""

    if config["mafft"]["quiet"] == "on":
        config["mafft"]["quiet"] = "--quiet"
    else:
        config["mafft"]["quiet"] = ""
fetchConfigParameters()

rule all:
    input:
        f"{os.path.splitext(sars2Ref)[0]}.msa",
        expand(
            f"{os.path.splitext(sars2Ref)[0]}.msa.{{ext}}",
            ext = ["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"]
        ),
        expand(
            f"/datater/wu/data/tmp/{{datadir}}_{{sample}}_{{kmerSize}}_allOligo.set",
            datadir = datadir,
            sample = sample,
            kmerSize = kmerSize
        ),
        expand(
            f"/datater/wu/data/tmp/{{datadir}}_{{sample}}_{{kmerSize}}_stritMatch.bowtie",
            datadir = datadir,
            sample = sample,
            kmerSize = kmerSize
        ),
        expand(
            f"/datater/wu/data/tmp/{{datadir}}_{{sample}}_{{kmerSize}}_notstritMatch.bowtie",
            datadir = datadir,
            sample = sample,
            kmerSize = kmerSize
        ),
        expand(
            f"/datater/wu/data/tmp/{{datadir}}_{{sample}}_{{kmerSize}}_statistic.out",
            datadir = datadir,
            sample = sample,
            kmerSize = kmerSize
        ),
        # expand(
        #     f"/datater/wu/data/tmp/{{datadir}}_{{sample}}_{{kmerSize}}_notstritMatch_statistic.out",
        #     datadir = datadir,
        #     sample = sample,
        #     kmerSize = kmerSize
        # )

rule test0:
    input:
        sars2Ref
    output:
        f"{os.path.splitext(sars2Ref)[0]}.msa"
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
        """
        mafft {params.algorithm} \
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
        {input} > {output}
        """

rule test1:
    input:
        f"{os.path.splitext(sars2Ref)[0]}.msa"
    output:
        multiext(
            f"{os.path.splitext(sars2Ref)[0]}.msa.", "1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"
        )
    params:
        out = f"{os.path.splitext(sars2Ref)[0]}.msa",
        thread = 12
    shell:
        """
        bowtie-build \
        --threads {params.thread} \
        {input[0]} {params.out}
        """

rule test2:
    input:
        f"/datater/wu/data/{{datadir}}/{{sample}}/{{kmerSize}}/allOligo.set",
    output:
        f"/datater/wu/data/tmp/{{datadir}}_{{sample}}_{{kmerSize}}_allOligo.set"
    run:
        import os

        outFile = open(output[0], "w")
        input1Open = open(input[0], "r")
        for li in input1Open.readlines()[1:]:
            li = li.rstrip("\n")
            ls = li.split()
            id = ls[0]
            oligo = ls[1]
            outFile.write(f">p{id}\n{oligo}\n")
        input1Open.close()
        outFile.close()

rule test3:
    input:
        f"/datater/wu/data/tmp/{{datadir}}_{{sample}}_{{kmerSize}}_allOligo.set",
        lambda wildcards: expand(
            f"{os.path.splitext(sars2Ref)[0]}.msa.{{ext}}",
            ext = ["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"]
        )
    output:
        f"/datater/wu/data/tmp/{{datadir}}_{{sample}}_{{kmerSize}}_stritMatch.bowtie"
    params:
        refDB = f"{os.path.splitext(sars2Ref)[0]}.msa",
        thread = 12
    shell:
        """
        bowtie \
        -x {params.refDB} -f {input[0]} \
        -v 0 -l 7 -a \
        --sam \
        -p {params.thread} \
        {output}
        """

rule test4:
    input:
        f"/datater/wu/data/tmp/{{datadir}}_{{sample}}_{{kmerSize}}_allOligo.set",
        lambda wildcards: expand(
            f"{os.path.splitext(sars2Ref)[0]}.msa.{{ext}}",
            ext = ["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"]
        )
    output:
        f"/datater/wu/data/tmp/{{datadir}}_{{sample}}_{{kmerSize}}_notstritMatch.bowtie"
    params:
        refDB = f"{os.path.splitext(sars2Ref)[0]}.msa",
        thread = 12
    shell:
        """
        bowtie \
        -x {params.refDB} -f {input[0]} \
        -n 2 -l 7 -a \
        --sam \
        -p {params.thread} \
        {output}
        """

rule test5:
    input:
        f"/datater/wu/data/tmp/{{datadir}}_{{sample}}_{{kmerSize}}_stritMatch.bowtie",
        f"/datater/wu/data/tmp/{{datadir}}_{{sample}}_{{kmerSize}}_notstritMatch.bowtie",
        f"{os.path.splitext(sars2Ref)[0]}.msa"
    output:
        f"/datater/wu/data/tmp/{{datadir}}_{{sample}}_{{kmerSize}}_statistic.out",
        # f"/datater/wu/data/tmp/{{datadir}}_{{sample}}_{{kmerSize}}_notstritMatch_statistic.out"
    params: 
        kmerSize = lambda wildcards: wildcards.kmerSize
    run:
        import os
        import pandas as pd
        import numpy as np
        from collections import defaultdict
        from Bio import SeqIO

        kmerSize = int(params.kmerSize.replace("kmer", ""))

        dict_strictMatch={}
        dict_notStrictMatch={}

        df_bowtieStrictMatch = pd.read_table(input[0], comment="@", sep="\t", names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "OPT1", "OPT2", "OPT3", "OPT4"])
        df_bowtieNotStrictMatch = pd.read_table(input[1], comment="@", sep="\t", names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "OPT1", "OPT2", "OPT3", "OPT4"])

        records = list(SeqIO.parse(input[2], "fasta"))
        genomeLen = len(records[0].seq)
        totalSequenceCount = len(records)
        chunks = [[i,i+500] for i in range(0, genomeLen, 500)]

        ### strict   
        df_bowtieStrictMatch = df_bowtieStrictMatch[df_bowtieStrictMatch["FLAG"]==0|16]
        df_bowtieStrictMatch = df_bowtieStrictMatch[["RNAME", "POS"]]
        df_bowtieStrictMatch = df_bowtieStrictMatch.rename(columns={"POS":"rstart"})
        df_bowtieStrictMatch["rend"] = df_bowtieStrictMatch["rstart"].apply(int) + kmerSize -1 
        df_bowtieStrictMatch["start"] = ""
        df_bowtieStrictMatch["end"] = ""

        criteria =[]
        for chunk in chunks:
            criteria.append(df_bowtieStrictMatch.rstart.between(chunk[0], chunk[1]))

        chunk_start = [x[0] for x in chunks]
        chunk_end = [x[1] for x in chunks]

        df_bowtieStrictMatch["start"] = np.select(criteria, chunk_start, 0)
        df_bowtieStrictMatch["end"] = np.select(criteria, chunk_end, 0)

        df_bowtieStrictMatch = df_bowtieStrictMatch.groupby(['start', 'end']).RNAME.nunique().reset_index(name='counts')
        # df_bowtieStrictMatch = df_bowtieStrictMatch.groupby(['start', 'end']).size().reset_index(name='counts')

        ### not strict
        df_bowtieNotStrictMatch = df_bowtieNotStrictMatch[df_bowtieNotStrictMatch["FLAG"]==0|16]
        df_bowtieNotStrictMatch = df_bowtieNotStrictMatch[["RNAME", "POS"]]
        df_bowtieNotStrictMatch = df_bowtieNotStrictMatch.rename(columns={"POS":"rstart"})
        df_bowtieNotStrictMatch["rend"] = df_bowtieNotStrictMatch["rstart"].apply(int) + kmerSize -1 
        df_bowtieNotStrictMatch["start"] = ""
        df_bowtieNotStrictMatch["end"] = ""

        criteria =[]
        for chunk in chunks:
            criteria.append(df_bowtieNotStrictMatch.rstart.between(chunk[0], chunk[1]))

        chunk_start = [x[0] for x in chunks]
        chunk_end = [x[1] for x in chunks]

        df_bowtieNotStrictMatch["start"] = np.select(criteria, chunk_start, 0)
        df_bowtieNotStrictMatch["end"] = np.select(criteria, chunk_end, 0)

        # df_bowtieNotStrictMatch = df_bowtieNotStrictMatch.groupby(['start', 'end']).size().reset_index(name='counts')
        df_bowtieNotStrictMatch = df_bowtieNotStrictMatch.groupby(['start', 'end']).RNAME.nunique().reset_index(name='counts')

        df_merge = pd.merge(df_bowtieStrictMatch, df_bowtieNotStrictMatch, on = ['start', 'end'], how = 'outer')
        df_merge = df_merge.fillna(0)
        df_merge["SeqCount"] = totalSequenceCount
        df_merge.to_csv(output[0], sep="\t", header=True, index=False)
        
