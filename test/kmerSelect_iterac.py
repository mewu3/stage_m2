#!/usr/bin/env python3.8
import sys
import os
import re

### global variables 
inputFile = sys.argv[1]
outputFile = sys.argv[2]
kmerCount = seqNumber = int(os.popen(f"egrep -c '^>' {inputFile}").read())
filename = os.path.splitext(inputFile)[0]
outputFile = open(outputFile, "w")

### functions 
def calculate_kmer(input):

    dskH5 = f"{filename}.h5"
    os.system(
        f"../lib/dsk/build/bin/dsk \
        -nb-cores 0 \
        -max-memory 5000 \
        -max-disk 0 \
        -out-compress 0 \
        -storage-type hdf5 \
        -kmer-size 13 \
        -abundance-min 2 \
        -abundance-max 2147483647 \
        -solidity-kind sum \
        -file {input} \
        -out {dskH5}"
    )

    dskh5txt = f"{filename}.kCount"
    os.system(
        f"../lib/dsk/build/bin/dsk2ascii \
        -nb-cores 0 \
        -file {dskH5} -out {dskh5txt}"
    )

    dskh5txtSorted = f"{filename}.kCountSorted"
    os.system(
        f"sort -s -n -k 2 -nr {dskh5txt} > {dskh5txtSorted}"
    )

    dict = {}

    file = open(dskh5txtSorted, "r") 
    for line in file.readlines()[0:2]:
        kmer = line.split()[0]
        kmerCount = line.split()[1]
        dict[kmer]=kmerCount
    file.close()

    return dict

def parse_fasta(input): 

    dict={}

    file = open(input, "r")

    for line in file: 
        line = line.rstrip("\n")
        if line.startswith(">"): 
            header = line 
            dict[header] = ""
        else: 
            dict[header] += line 

    file.close()

    return(dict)

while kmerCount > seqNumber*0.1: 
   
    kmerDict = calculate_kmer(inputFile)
    kmerList = [key for key in kmerDict if kmerDict]     
    fastaDict = parse_fasta(inputFile)

    if kmerList: 
        for kmer in kmerList: 
            boolean = [bool(re.search(kmer, fastaDict[key], re.I)) for key in fastaDict]
            if True in boolean: 
                kmerCount = int(kmerDict[kmer])
                outputFile.write(f"{kmer}\t{kmerDict[kmer]}\n")
                intermediate = open(f"{filename}.inter", "w")
                for key in fastaDict: 
                    seq = fastaDict[key]
                    if re.search(kmer, seq, re.I) is None: 
                        intermediate.write(f"{key}\n{seq}\n")
                intermediate.close()
                break 
            
    inputFile = f"{filename}.inter"

    # ~ boolean = [bool(re.search(kmerList[0], fastaDict[key], re.I)) for key in fastaDict]

    # ~ if True in boolean: 
        # ~ outputFile.write(f"{kmerList[0]}\t{kmerDict[kmerList[0]]}\n")
        # ~ intermediate = open(f"{filename}.inter", "w")
        # ~ for key in fastaDict: 
            # ~ seq = fastaDict[key]
            # ~ if re.search(kmerList[0], seq, re.I) is None: 
                # ~ intermediate.write(f"{key}\n{seq}\n")
        # ~ intermediate.close()
    # ~ else: 
        # ~ outputFile.write(f"{kmerList[1]}\t{kmerDict[kmerList[1]]}\n")
        # ~ intermediate = open(f"{filename}.inter", "w")
        # ~ for key in fastaDict: 
            # ~ seq = fastaDict[key]
            # ~ if re.search(kmerList[1], seq, re.I) is None: 
                # ~ intermediate.write(f"{key}\n{seq}\n")
        # ~ intermediate.close()
    
    # ~ inputFile = f"{filename}.inter"

outputFile.close()
