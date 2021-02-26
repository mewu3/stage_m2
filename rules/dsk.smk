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
        histoMax = config["dsk"]["histo-max"],# with os.scandir("test/splitFiles") as it:
#     for entry in it:
#         if not entry.name.startswith('.') and entry.is_file():
#             print(entry.name)
        histo2D = config["dsk"]["histo2D"],
        histo = config["dsk"]["histo"]
    script:
        "scripts/dsk.py"
