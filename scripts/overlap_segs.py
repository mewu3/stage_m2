#!/usr/bin/env python3.8
import sys
import os 
import numpy as np
import pandas as pd
from pyfaidx import Fasta


### input arguments
inputFile = snakemake.input.fasta
#inputFile = sys.argv[1]
#step = sys.argv[2]
step = 500
#overlap = sys.argv[3]
overlap = 250
#region = sys.argv[4]
outputFile = snakemake.output.dir
os.makedirs(outputFile)
#outputFile = sys.argv[2]

### global variables
fasta = Fasta(inputFile)
seqLength = len(fasta[0])

### split into overlapping segments
remainder = seqLength % (step-overlap)
chunkNumber = int(seqLength / (step-overlap))
print("The sequences are splitted into "+str(chunkNumber)+" chunks, and there are "+str(remainder)+" bp left.")


chunks = [[i,i+step] for i in range(0, seqLength, step-overlap)]
chunks[-1][1] = len(fasta[0])


for chunk in chunks: 
    f = open(outputFile+"/"+"forwardChunk"+str(chunk[0])+"-"+str(chunk[1])+".fasta", "w")
    for id in fasta.keys() : 
        segment = fasta[id][chunk[0]:chunk[1]]
        forward = str(segment[:50])
        reverse = str(segment[-50:])
        f.write(">" + fasta[id].long_name + " |" + str(chunk[0]) + "-" + str(chunk[1])+"\n")
        f.write(forward+"\n")
    f.close()

