#!/usr/bin/env python3.8
import sys
from pyfaidx import Fasta
import os

inputFile = sys.argv[1]
# os.makedirs(output[0])
step = int(sys.argv[2])
overlap = int(sys.argv[3])

fasta = Fasta(inputFile)
seqLength = len(fasta[0])

remainder = seqLength % (step-overlap)
chunkNumber = int(seqLength / (step-overlap))
print("The sequences are splitted into "+str(chunkNumber)+" chunks, and there are "+str(remainder)+" bp left.")

if remainder <= step/2: # primux fasta_tile_overlap.pl
    newStep = int(remainder/chunkNumber) + 1 + step
    print("Changing step size from {} to {} so there will be no remainder.".format(step, newStep))
    step = newStep


chunks = [[i,i+step] for i in range(0, seqLength, step-overlap)]
chunkNumber=len(chunks)

for chunk in chunks:

    if chunk[-1] > seqLength:
        chunk[-1] = seqLength

    seg = str(chunk[0])+"-"+str(chunk[1])
    print(seg)

    # f1 = open(output[0] + f"/forward{seg}.fasta", "w")
    # f2 = open(output[0] + f"/reverse{seg}.fasta", "w")

    for id in fasta.keys() :
        segment = fasta[id][chunk[0]:chunk[1]]
        forward = str(segment[:50])
        reverse = str(segment[-50:])
        # f1.write("> {}|{}-{}\n{}\n".format(fasta[id].long_name, str(chunk[0]), str(chunk[1]), forward))
        # f2.write("> {}|{}-{}\n{}\n".format(fasta[id].long_name, str(chunk[0]), str(chunk[1]), reverse))

    # f1.close()
    # f2.close()

print(f"Final chunk number: {chunkNumber}")
