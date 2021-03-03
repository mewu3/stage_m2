rule dsk: # Kmer counting ######################################################
    input:
        datadir + "/splitFiles/{prefix}.fasta"
    output:
        temp(datadir + "/kmerCounting/dsk/{prefix}.h5")
    params:
        nbCores = config["dsk"]["nb-cores"],
        maxMemory = config["dsk"]["max-memory"],
        maxDisk = config["dsk"]["max-disk"],
        outCompress = config["dsk"]["out-compress"],
        storage = config["dsk"]["storage-type"],
        kmerSize = config["dsk"]["kmer-size"],
        abundanceMin = config["dsk"]["abundance-min"],
        abundanceMax = config["dsk"]["abundance-max"],
        solidityKind = config["dsk"]["solidity-kind"]
    shell:
        "lib/dsk/build/bin/dsk \
        -nb-cores {params.nbCores} \
        -max-memory {params.maxMemory} \
        -max-disk {params.maxDisk} \
        -out-compress {params.outCompress} \
        -storage-type {params.storage} \
        -kmer-size {params.kmerSize} \
        -abundance-min {params.abundanceMin} \
        -abundance-max {params.abundanceMax} \
        -solidity-kind {params.solidityKind} \
        -file {input} -out {output}"

rule dskOutput:
    input:
        datadir + "/kmerCounting/dsk/{prefix}.h5"
    output:
        datadir + "/kmerCounting/dsk/{prefix}.txt"
    params:
        nbCores = config["dsk2ascii"]["nb-cores"],
    shell:
        "lib/dsk/build/bin/dsk2ascii \
        -nb-cores {params.nbCores} \
        -file {input} -out {output}"

rule dskOutputSort:
    input:
        datadir + "/kmerCounting/dsk/{prefix}.txt"
    output:
        datadir + "/kmerCounting/dsk/{prefix}.sort"
    run:
        #!/usr/bin/env python3.8
        with open(input[0], "r") as input:
            rows = input.readlines()
            sorted_rows = sorted(rows, key = lambda x: int(x.split()[1]), reverse=True)
            with open(output[0], "w") as output:
                for row in sorted_rows:
                        output.write(row)
