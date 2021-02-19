#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
from pyfasta import Fasta

### input arguments
input = sys.argv[1]
# step = sys.argv[2]
# overStep = sys.argv[3]
# region = sys.argv[4]

input = Fasta(input)
seqNum_count = len(input)
