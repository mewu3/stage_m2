rule removeDuplicateSeq:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        f"{dataDir}/{{sample}}/{{sample}}.uniq"
    shell:
        "lib/cd-hit-v4.8.1-2019-0228/cd-hit-auxtools/cd-hit-dup -i {input} -o {output}"

rule countAllKmer1:
    input:
        f"{dataDir}/{{sample}}/{{sample}}.uniq"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.jf"
    params:
        kmerSize = config["jellyfish"]["kmer-size"],
        hash = config["jellyfish"]["hash-size"],
        thread = config["jellyfish"]["thread"]
    shell:
        """
        jellyfish count \
        -m {params.kmerSize} \
        -s {params.hash} \
        -t {params.thread} \
        -o {output} {input}
        """

rule countAllKmer2:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.jf"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.txt"
    params:
        lowerCount = config["jellyfish"]["lower-count"]
    shell:
        "jellyfish dump -c -t -L {params.lowerCount} {input} > {output}"

rule sortAllKmer:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.txt"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.txt"
    shell:
        "cat {input} | sort -n -k2 -r > {output}"
