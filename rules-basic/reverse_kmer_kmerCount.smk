rule countKmer1:
    input:
        f"{dataDir}/{{sample}}/splitFiles/reverse{{seg}}.fasta"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/reverse{{seg}}.jf"
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
        -o {output} {input} \
        """

rule countKmer2:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/reverse{{seg}}.jf"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/reverse{{seg}}.Kmercount.txt"
    shell:
        """
        jellyfish dump -c {input} > {output}
        """

rule calculateKmer_sort:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/reverse{{seg}}.Kmercount.txt"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/reverse{{seg}}.Kmercount.sorted.txt"
    shell:
        "cat {input} | sort -n -k2 -r > {output}"
