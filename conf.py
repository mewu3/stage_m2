#!/usr/bin/env python
import pandas as pd
import numpy as np
import os



sample = pd.read_table("samples.tsv").set_index("sample", drop=False)
print(sample)
