#!/usr/bin/env python
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
lenRange = range(0,len(input[0]))

### split into overlapping segments
chunks = [lenRange[i:i+step] for i in range(0, len(input[0]), step-overlap)]

print(chunks[0])
print(input[0][chunks[0]])
