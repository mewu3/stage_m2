#!/usr/bin/env python3.8
import sys
import os

inputfile = sys.argv[1]
# outputfile = open(sys.argv[2], "w")
# dataDir = snakemake.params.dataDir
dataDir = "/datater/wu/data"
sample = "enterovirus"
seg = "0-504"
# sample = snakemake.wildcards.sample
# seg = snakemake.wildcards.seg



def calculate_kmer(inputfile):
    dskH5 = f"{dataDir}/{sample}/dsk/{seg}.h5"
    os.system(
        f"lib/dsk/build/bin/dsk \
        -nb-cores 0 \
        -max-memory 5000 \
        -max-disk 0 \
        -out-compress 0 \
        -storage-type hdf5 \
        -kmer-size 13 \
        -abundance-min 2 \
        -abundance-max 2147483647 \
        -solidity-kind sum \
        -file {inputfile} \
        -out {dskH5}"
    )

    dskh5txt = f"{dataDir}/{sample}/dsk/{seg}.txt"
    os.system(
        f"lib/dsk/build/bin/dsk2ascii \
        -nb-cores 0 \
        -file {dskH5} -out {dskh5txt}"
    )

    dskh5txtSorted = f"{dataDir}/{sample}/dsk/{seg}.sorted"
    os.system(
        f"sort -s -n -k 2 -nr {dskh5txt} > {dskh5txtSorted}"
    )

    kmer1st = os.system(f"head -n 1 {dskh5txtSorted}")
    return kmer1st

calculate_kmer(inputfile)
