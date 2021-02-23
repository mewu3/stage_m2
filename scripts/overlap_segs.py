#!/usr/bin/env python3.8
import sys
import numpy as np
import pandas as pd
from pyfaidx import Fasta

### input arguments
input = sys.argv[1]
# step = sys.argv[2]
step = 500
# overlap = sys.argv[3]
overlap = 250
# region = sys.argv[4]

### global variables
input = Fasta(input)
lsRange = list(range(len(input[0])))
ran = range(len(input[0]))

### split into overlapping segments
chunks = [[i,i+step] for i in range(0, len(input[0]), step-overlap)]

# for id in input.keys(): 
#     allFoward = []
#     allReverse = []
#     for chunk in chunks:
#         segment = input[id][chunk[0]:chunk[1]]
#         forward = str(segment[:50])
#         reverse = str(segment[-50:])
#         allFoward.append(forward)
#     print(allFoward)

for chunk in chunks: 
    for id in input.keys() : 
        segment = input[id][chunk[0]:chunk[1]]
        forward = str(segment[:50])
        reverse = str(segment[-50:])
        print(forward, end="\n")
    print()