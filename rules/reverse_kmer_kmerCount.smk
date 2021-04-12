rule calculateKmer_jellyfishCount:
    input:
        f"{dataDir}/{{sample}}/splitFiles/reverse{{seg}}.fasta"
    output:
        f"{dataDir}/{{sample}}/kmerCounting{kmerSize}/reverse{{seg}}.jf"
    params:
        kmerSize = config["jellyfish-count"]["kmer-size"],
        hash = config["jellyfish-count"]["hash-size"],
        thread = config["jellyfish-count"]["thread"]
    shell:
        "jellyfish count \
        -m {params.kmerSize} \
        -s {params.hash} \
        -t {params.thread} \
        -o {output} {input}"

rule calculateKmer_jellyfishDump:
    input:
        f"{dataDir}/{{sample}}/kmerCounting{kmerSize}/reverse{{seg}}.jf"
    output:
        f"{dataDir}/{{sample}}/kmerCounting{kmerSize}/reverse{{seg}}.kCount"
    shell:
        "jellyfish dump -c {input} > {output}"

rule calculateKmer_sort:
    input:
        f"{dataDir}/{{sample}}/kmerCounting{kmerSize}/reverse{{seg}}.kCount"
    output:
        f"{dataDir}/{{sample}}/kmerCounting{kmerSize}/reverse{{seg}}.kCountSorted"
    shell:
        "cat {input} | sort -n -k2 -r > {output}"
