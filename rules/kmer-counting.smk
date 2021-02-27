rule dsk:
    input:
        dir = "test/splitFiles/"
    output:
        dir = directory("test/kmerCounting/dsk")
    params:
        nbCores = config["dsk"]["nb-cores"],
        maxMemory = config["dsk"]["max-memory"],
        maxDisk = config["dsk"]["max-disk"],
        outCompress = config["dsk"]["out-compress"],
        storage = config["dsk"]["storage-type"],
        verbose = config["dsk"]["verbose"],
        kmerSize = config["dsk"]["kmer-size"],
        abundanceMin = config["dsk"]["abundance-min"],
        abundanceMax = config["dsk"]["abundance-max"],
        abundanceMinThreshold = config["dsk"]["abundance-min-threshold"],
        solidityKind = config["dsk"]["solidity-kind"],
        solidityCustom = config["dsk"]["solidity-custom"],
        solideKmerOut = config["dsk"]["solid-kmers-out"],
        histoMax = config["dsk"]["histo-max"],
        histo2D = config["dsk"]["histo2D"],
        histo = config["dsk"]["histo"]
    script:
        "scripts/dsk.py"

rule dskOutput:
    input:
        dir = "test/kmerCounting/dsk"
    output:
        dir = directory("test/kmerCounting/dsk_output")
    script:
        "scripts/dsk_output.py"

rule kmc:
    input:
        dir = "test/splitFiles"
    output:
        dir = directory("test/kmerCounting/kmc3")
    # log:
    #     "test/kmerCounting/kmc3/log.txt"
    conda:
        "envs/kmc3.yaml"
    params:
        kmerSize = config["kmc3"]["kmer-size"]
    script:
        "scripts/kmc3.py"

rule kmcOutput:
    input:
        dir = "test/kmerCounting/kmc3"
    output:
        dir = directory("test/kmerCounting/kmc3_output")
    script:
        "scripts/kmc3_output.py"
    script:
        "scripts/dsk.py"
