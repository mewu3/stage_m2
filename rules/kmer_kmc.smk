rule kmc: # Kmer counting ######################################################
    input:
        datadir + "/splitFiles/{prefix}.fasta"
    output:
        datadir + "/kmerCounting/kmc/{prefix}.kmc_pre",
        datadir + "/kmerCounting/kmc/{prefix}.kmc_suf"
    conda:
        "envs/kmc3.yaml"
    params:
        prefix = datadir + "/kmerCounting/kmc/{prefix}",
        outdir = datadir + "/kmerCounting/kmc",
        kmerSize = config["kmc3"]["kmer-size"],
        inputFormat = config["kmc3"]["input-format"],
        ram = config["kmc3"]["ram"],
        signature = config["kmc3"]["signature"],
        minKmer = config["kmc3"]["min-kmer-number"],
        maxKmer = config["kmc3"]["max-kmer-number"],
        maxKmerExclu = config["kmc3"]["max-exclu-kmer"]
    shell:
        "kmc -k{params.kmerSize} \
        -m{params.ram} \
        -cs{params.maxKmer} \
        {params.inputFormat} {input} \
        {params.prefix} {params.outdir}"

rule kmcOutput:
    input:
        datadir + "/kmerCounting/kmc/{prefix}.kmc_pre",
        datadir + "/kmerCounting/kmc/{prefix}.kmc_suf"
    output:
        datadir + "/kmerCounting/kmc/{prefix}.txt"
    params:
        prefix = datadir + "/kmerCounting/kmc/{prefix}"
    shell:
        "kmc_dump {params.prefix} {output}"

rule kmcOutputSort:
    input:
        datadir + "/kmerCounting/kmc/{prefix}.txt"
    output:
        datadir + "/kmerCounting/kmc/{prefix}.sort"
    run:
        #!/usr/bin/env python3.8
        with open(input[0], "r") as input:
            rows = input.readlines()
            sorted_rows = sorted(rows, key = lambda x: int(x.split()[1]), reverse=True)
            with open(output[0], "w") as output:
                for row in sorted_rows:
                        output.write(row)
