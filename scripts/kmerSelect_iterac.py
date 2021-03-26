#!/usr/bin/env python3.8
import sys
import os
from Bio import SeqIO
import re
from collections import Counter

inputfile = sys.argv[1]
outputfile = open(sys.argv[2], "w")
dataDir = "/datater/wu/data"
sample = "enterovirus"
seg = "0-504"
seqNum = int(os.popen(f"grep '^>' -c {inputfile}").read())
# dataDir = snakemake.params.dataDir
# sample = snakemake.wildcards.sample
# seg = snakemake.wildcards.seg

def calculate_kmer(inputfile):
    dskH5 = f"{dataDir}/{sample}/dsk/{seg}.inter.h5"
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

    dskh5txt = f"{dataDir}/{sample}/dsk/{seg}.inter.txt"
    os.system(
        f"lib/dsk/build/bin/dsk2ascii \
        -nb-cores 0 \
        -file {dskH5} -out {dskh5txt}"
    )

    dskh5txtSorted = f"{dataDir}/{sample}/dsk/{seg}.inter.sorted"
    os.system(
        f"sort -s -n -k 2 -nr {dskh5txt} > {dskh5txtSorted}"
    )

    kmerDict = {}
    with open(dskh5txtSorted, "r") as file:
        for line in file.readlines()[0:3]:
            kmer = line.split()[0]
            kmerCount = line.split()[1]
            kmerDict[kmer]=kmerCount
    return kmerDict

n=0
kmerCount = seqNum



# for n in range(0,3):
# # # while int(kmerCount) > int(seqNum*0.1):
#     kmerList=[]
#     kmerDict = calculate_kmer(inputfile)
#     kmerList = [key for key in kmerDict]
#     record_dict = SeqIO.to_dict(SeqIO.parse(inputfile, "fasta"))
#
#     boolean = [bool(re.search(kmerList[0], str(record_dict[record].seq), re.IGNORECASE)) for record in record_dict]
#
#     for kmer in kmerList:
#
#         if True in boolean:
#             intermediate = f"{dataDir}/{sample}/dsk/{seg}.inter{n}.truncated"
#             intermediate = open(intermediate, "w")
#             for record in record_dict:
#                 seq = str(record_dict[record].seq)
#                 if re.search(kmer, seq, re.IGNORECASE) is None:
#                     outputfile.write(f"{kmer}\t{kmerDict[kmer]}\n")
#                 else:
#                     intermediate.write(record_dict[record].format("fasta"))
#                     outputfile.write(f"{kmer}\t{kmerDict[kmer]}\n")
#             intermediate.close()
#     inputfile = intermediate
#     n += 1

for n in range(0,3):
    kmerList=[]
    kmerDict = calculate_kmer(inputfile)
    kmerList = [key for key in kmerDict]
    record_dict = SeqIO.to_dict(SeqIO.parse(inputfile, "fasta"))
    intermediate = open(f"{dataDir}/{sample}/dsk/{seg}.inter{n}.truncated", "w")

    for kmer in kmerList[0:3]:
        for record in record_dict:
            seq = str(record_dict[record].seq)
            boolean = [bool(re.search(kmer, seq, re.IGNORECASE)) for record in record_dict]
            if True in boolean:
                if re.search(kmer, seq, re.IGNORECASE) is None:
                    outputfile.write(f"{kmer}\t{kmerDict[kmer]}\n")
                    intermediate.write(record_dict[record].format("fasta"))
    intermediate.close()
    inputfile = f"{dataDir}/{sample}/dsk/{seg}.inter{n}.truncated"
    n+=1

    outputfile.close()
