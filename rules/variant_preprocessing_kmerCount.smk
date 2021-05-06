rule removeDuplicateSeq:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        f"{dataDir}/{{sample}}/{{sample}}.uniq"
    shell:
        "lib/cd-hit-v4.8.1-2019-0228/cd-hit-auxtools/cd-hit-dup -i {input} -o {output}"

# rule countAllKmer1:
#     input:
#         f"{dataDir}/{{sample}}/{{sample}}.uniq"
#     output:
#         f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.int"
#     params:
#         kmerSize = config["kmerSize"],
#         hash = config["jellyfish"]["hash-size"],
#         thread = config["jellyfish"]["thread"]
#     shell:
#         """
#         jellyfish count \
#         -m {params.kmerSize} \
#         -s {params.hash} \
#         -t {params.thread} \
#         -o {output} {input}
#         """
#
# rule countAllKmer2:
#     input:
#         f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.int"
#     output:
#         f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.txt"
#     shell:
#         "jellyfish dump -c {input} > {output}"
#
# rule sortAllKmer:
#     input:
#         f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.txt"
#     output:
#         f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.txt"
#     shell:
#         "cat {input} | sort -n -k2 -r > {output}"

rule countAllKmer1:
    input:
        lambda wildcards: config["samples"][wildcards.sample] if config["curated"] else f"{dataDir}/{{sample}}/{{sample}}.uniq"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.int.kmc_pre",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.int.kmc_suf",
    params:
        kmerSize = lambda wildcards: config["kmerSize"][wildcards.kmerSize],
        memory = config["kmc3"]["memory"],
        maxCount = config["kmc3"]["maxCount"],
        thread = config["thread"],
        output = f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.int",
        tmpdir = f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate"
    shell:
        """kmc -k{params.kmerSize} \
        -m{params.memory} -b \
        -cs{params.maxCount} \
        -r{params.thread} \
        -fm {input[0]} {params.output} {params.tmpdir}"""

rule countAllKmer2:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.int.kmc_pre",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.int.kmc_suf",
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.txt"
    params:
        output = f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.int"
    shell:
        "kmc_dump {params.output} {output[0]}"

rule sortAllKmer:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.txt"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.txt"
    shell:
        "cat {input} | sort -n -k2 -r > {output}"
