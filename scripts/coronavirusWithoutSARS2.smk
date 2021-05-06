# !/usr/bin/env python3.8
import os

configfile: "config.yaml"
samples = list(config["samples"].keys())
file_refSeq = config["file_refSeq"]
file_aceIDtaxID = config["file_aceIDtaxID"]
file_taxIDLineage = config["file_taxIDLineage"]
dataDir = config["dataDir"]
kmerSize = list(config["kmerSize"].keys())

rule all:
    input:
        "/datater/wu/data/tmp/sars2Complet.uniq.fasta",
        "/datater/wu/data/tmp/sars2Complet.uniq.cls95.fasta",
        # expand(
        #     f"/datater/wu/data/tmp/sars2Complet.uniq.cls95.fasta.DB.{{ext}}",
        #     ext = ["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"]
        # ),
        # expand(
        #     f"/datater/wu/data/tmp/{{sample}}/{{kmerSize}}/stritMatch.bowtie",
        #     sample = samples,
        #     kmerSize = kmerSize,
        # )

rule removeDuplicateSeq:
    input:
        "/datater/wu/data/tmp/sars2Complet.fasta"
    output:
        "/datater/wu/data/tmp/sars2Complet.uniq.fasta"
    conda:
        "envs/seqkit.yaml"
    shell:
        "lib/cdhit/cd-hit-auxtools/cd-hit-dup -i {input} -o {output}"

rule clustering:
    input:
        "/datater/wu/data/tmp/sars2Complet.uniq.fasta"
    output:
        "/datater/wu/data/tmp/sars2Complet.uniq.cls95.fasta"
    params:
        identity = config["cd-hit"]["identity"],
        threads = config["thread"],
        memory = config["cd-hit"]["memory"]
    shell:
        """
        lib/cdhit/psi-cd-hit//psi-cd-hit.pl -i {input[0]} -o {output[0]} \
        -c 0.95 -G 1 -g 1 -prog megablast -s "-F F -e 0.000001 -b 100000 -v 100000" -exec local -para 8 -blp 4
        """

rule test1:
    input:
        "/datater/wu/data/tmp/sars2Complet.uniq.cls95.fasta"
    output:
        multiext(
            f"/datater/wu/data/tmp/sars2Complet.uniq.cls95.fasta.DB.", "1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"
        )
    params:
        out = "/datater/wu/data/tmp/sars2Complet.uniq.cls95.fasta.DB",
        thread = 12
    shell:
        """
        bowtie-build \
        --threads {params.thread} \
        {input[0]} {params.out}
        """

rule test2:
    input:
        f"/datater/wu/data/MSSPE-variant/{{sample}}/{{kmerSize}}/allOligo.set",
        lambda wildcards: expand(
            f"/datater/wu/data/tmp/sars2Complet.uniq.cls95.fasta.DB.{{ext}}",
            ext = ["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"]
        )
    output:
        "/datater/wu/data/tmp/{{sample}}/{{kmerSize}}/stritMatch.bowtie"
    params:
        refDB = "/datater/wu/data/tmp/sars2Complet.uniq.cls95.fasta.DB",
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
