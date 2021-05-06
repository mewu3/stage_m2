rule countKmer1:
    input:
        f"{dataDir}/{{sample}}/splitFiles/reverse{{seg}}.fasta"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/reverse{{seg}}.int.kmc_pre",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/reverse{{seg}}.int.kmc_suf"
    params:
        kmerSize = lambda wildcards: config["kmerSize"][wildcards.kmerSize],
        memory = config["kmc3"]["memory"],
        maxCount = config["kmc3"]["maxCount"],
        thread = config["thread"],
        output = f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/reverse{{seg}}.int",
        tmpdir = f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate"
    shell:
        """kmc -k{params.kmerSize} \
        -m{params.memory} -b \
        -cs{params.maxCount} \
        -r{params.thread} \
        -fm {input[0]} \
        {params.output} {params.tmpdir}"""

rule countKmer2:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/reverse{{seg}}.int.kmc_pre",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/reverse{{seg}}.int.kmc_suf"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/reverse{{seg}}.Kmercount.txt"
    params:
        output = f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/reverse{{seg}}.int"
    shell:
        "kmc_dump {params.output} {output[0]}"

rule calculateKmer_sort:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/reverse{{seg}}.Kmercount.txt"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/reverse{{seg}}.Kmercount.sorted.txt"
    shell:
        "cat {input} | sort -n -k2 -r > {output}"
