#!/usr/bin/env python3.8

import os
import sys

os.makedirs(snakemake.output.dir)

with os.scandir(snakemake.input[0]) as it:
    for entry in it:
        if not entry.name.startswith('.') and entry.is_file():
            filename = os.path.splitext(entry.name)[0]
            os.system("lib/dsk/build/bin/dsk2ascii -file "+entry.path+" -out "+snakemake.output.dir+"/"+filename+".txt")
