#/usr/bin/env python3.8
import os
import sys

os.makedirs(snakemake.output.dir)

with os.scandir(snakemake.input[0]) as it:
    for entry in it:
        if not entry.name.startswith('.') and entry.is_file():
            # shell("lib/dsk/build/bin/dsk -nb-cores {params.nbCores} -max-memory {params.maxMemory} -max-disk {params.maxDisk} -out-compress {params.outCompress} -storage-type {params.storage} -verbose {params.verbose} -kmer-size {params.kmerSize} -abundance-min {params.abundanceMin} -abundance-max {params.abundanceMax} -abundance-min-threshold {params.abundanceMinThreshold} -solidity-kind {params.solidityKind} -solidity-custom {params.solidityCustom} -solid-kmers-out {params.solideKmerOut} -histo-max {params.histoMax} -histo2D {params.histo2D} -histo {params.histo} -file " + entry.path + "-out-dir " + output.dir)
            os.system("lib/dsk/build/bin/dsk -file " + entry.path + " -out-dir " + snakemake.output.dir)
