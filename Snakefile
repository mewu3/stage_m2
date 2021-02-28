#!/usr/bin/env python3.8

configfile: "config.yaml"

sample = list(config["samples"].keys())
outputdir = config["outputdir"]

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

rule all: # run all rules ######################################################
    input:
        expand(
            outputdir+"/duplicate_removed/{sample}.uniq.fasta",
            sample = sample
        ),
        expand(
            outputdir + "/msa/{sample}.msa.fasta",
               sample = sample
        ),
        dynamic(
            expand(
                outputdir + "/splitFiles/{sample}.forward.{seg}.fasta",
                sample = sample,
                seg="{seg}"
            )
        ),
        dynamic(
            expand(
                outputdir + "/splitFiles/{sample}.reverse.{seg}.fasta",
                sample = sample,
                seg="{seg}"
            )
        )


rule seqkit: # remove duplicate sequences ######################################
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        # "{}/duplicate_removed/{sample}.uniq.fasta".format(outputdir)
        outputdir + "/duplicate_removed/{sample}.uniq.fasta"
    conda:
        "envs/seqkit.yaml"
    log:
        # "{}/duplicate_removed/seqkit.{sample}.log".format(outputdir)
        outputdir + "/duplicate_removed/{sample}.uniq.fasta"
    shell:
        "cat {input} | seqkit rmdup -s -o {output}"

rule mafft: # multiple sequences alignment #####################################
    input:
        outputdir + "/duplicate_removed/{sample}.uniq.fasta"
    output:
        outputdir + "/msa/{sample}.msa.fasta"
    log:
        outputdir + "/msa/mafft.{sample}.log"
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
        "mafft \
        --thread {params.threads} \
        {params.algorithm} \
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

rule split_overlap_chunks: # split MSA fasta file into seperates files #########
    input:
        outputdir + "/msa/{sample}.msa.fasta"
    output:
        dynamic(outputdir + "/splitFiles/{sample}.forward.{seg}.fasta"),
        dynamic(outputdir + "/splitFiles/{sample}.reverse.{seg}.fasta")
    params:
        step = config["step"],
        overlap = config["overlap"]
    run:
        #!/usr/bin/env python3.8
        import sys
        import os
        import numpy as np
        import pandas as pd
        from pyfaidx import Fasta

        inputFile = input[0]
        step = params.step
        overlap = params.overlap

        fasta = Fasta(inputFile)
        seqLength = len(fasta[0])

        remainder = seqLength % (step-overlap)
        chunkNumber = int(seqLength / (step-overlap))
        print("The sequences are splitted into "+str(chunkNumber)+" chunks, and there are "+str(remainder)+" bp left.")

        if remainder <= step/2: # primux fasta_tile_overlap.pl
            newStep = int(remainder/chunkNumber) + 1 + step
            print("Changing step size from {} to {} so there will be no remainder.".format(step, newStep))

        chunks = [[i,i+newStep] for i in range(0, seqLength, newStep-overlap)]
        chunks[-1][1] = len(fasta[0])

        for chunk in chunks:

            seg = str(chunk[0])+"-"+str(chunk[1])

            f1 = open(outputdir + "/splitFiles/"+wildcards.sample+".forward.{}.fasta".format(seg), "w")
            f2 = open(outputdir + "/splitFiles/"+wildcards.sample+".reverse.{}.fasta".format(seg), "w")

            for id in fasta.keys() :
                segment = fasta[id][chunk[0]:chunk[1]]
                forward = str(segment[:50])
                reverse = str(segment[-50:])
                f1.write("> {} |{}-{} \n {} \n".format(fasta[id].long_name, str(chunk[0]), str(chunk[1]), forward))
                f2.write("> {} |{}-{} \n {} \n".format(fasta[id].long_name, str(chunk[0]), str(chunk[1]), reverse))

            f1.close()
            f2.close()



# rule dsk: # Kmer counting ######################################################
#     input:
#         dir = "test/splitFiles/"
#     output:
#         dir = directory("test/kmerCounting/dsk")
#     params:
#         nbCores = config["dsk"]["nb-cores"],
#         maxMemory = config["dsk"]["max-memory"],
#         maxDisk = config["dsk"]["max-disk"],
#         outCompress = config["dsk"]["out-compress"],
#         storage = config["dsk"]["storage-type"],
#         verbose = config["dsk"]["verbose"],
#         kmerSize = config["dsk"]["kmer-size"],
#         abundanceMin = config["dsk"]["abundance-min"],
#         abundanceMax = config["dsk"]["abundance-max"],
#         abundanceMinThreshold = config["dsk"]["abundance-min-threshold"],
#         solidityKind = config["dsk"]["solidity-kind"],
#         solidityCustom = config["dsk"]["solidity-custom"],
#         solideKmerOut = config["dsk"]["solid-kmers-out"],
#         histoMax = config["dsk"]["histo-max"],
#         histo2D = config["dsk"]["histo2D"],
#         histo = config["dsk"]["histo"]
#     script:
#         "scripts/dsk.py"
#
# rule dskOutput:
#     input:
#         dir = "test/kmerCounting/dsk"
#     output:
#         dir = directory("test/kmerCounting/dsk_output")
#     script:
#         "scripts/dsk_output.py"
#
# rule kmc: # Kmer counting ######################################################
#     input:
#         dir = "test/splitFiles"
#     output:
#         dir = directory("test/kmerCounting/kmc3")
#     # log:
#     #     "test/kmerCounting/kmc3/log.txt"
#     conda:
#         "envs/kmc3.yaml"
#     params:
#         kmerSize = config["kmc3"]["kmer-size"]
#     script:
#         "scripts/kmc3.py"
#
# rule kmcOutput:
#     input:
#         dir = "test/kmerCounting/kmc3"
#     output:
#         dir = directory("test/kmerCounting/kmc3_output")
#     script:
#         "scripts/kmc3_output.py"
#     script:
#         "scripts/dsk.py"
