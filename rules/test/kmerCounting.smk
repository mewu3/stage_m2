rule kmerCounting1:
    input:
        input_kmerCounting1
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
        """
        kmc -k{params.kmerSize} \
        -m{params.memory} -b \
        -cs{params.maxCount} \
        -r{params.thread} \
        -fm {input[0]} {params.output} {params.tmpdir}
        """

rule kmerCounting2:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.int.kmc_pre",
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.int.kmc_suf",
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.txt"
    params:
        output = f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.int"
    shell:
        "kmc_dump {params.output} {output[0]}"

rule kmerCounting3:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.txt"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.txt"
    shell:
        "cat {input} | sort -n -k2 -r > {output}"