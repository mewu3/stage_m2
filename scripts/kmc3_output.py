#/usr/bin/env python3.8

import os

os.makedirs(snakemake.output.dir)

with os.scandir(snakemake.input[0]) as it:
    for entry in it:
        if not entry.name.startswith(".") and entry.is_file():
            filename = os.path.splitext(entry.name)[0]
            noextention = entry.path.split(".")[0]
            # print(noextention)
            os.system("kmc_dump -ci2 -cx1000 "+noextention+" "+snakemake.output.dir+"/"+filename+".txt")