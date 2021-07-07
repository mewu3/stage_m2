import sys, os
from Bio import SeqIO
from pathlib import Path
import pandas as pd
import glob

input = "/datater/wu/data/test_clustering/enterovirus_all.cluster.clstr"

# with open(input, "r") as f:
#     for l in f:
#         l = l.rstrip("\n")
#         if l.startswith(">"):
#             cluster = l
#             cluster = l.lstrip(">").replace(" ", "_")
#             dir = f"/datater/wu/data/test_clustering/enterovirus/{cluster}/"
#             Path(dir).mkdir(parents=True, exist_ok=True)
#         else:
#             l=l.split()
#             aceID=l[2].lstrip(">").rstrip("...")
#             os.system(f"mv /datater/wu/data/test_clustering/enterovirus/{aceID}.fasta {dir}")

# for clus in $(ls ./ ); do realpath $clus/*fasta >> $clus/${clus}.txt ; done

# for fasta in $(ls Cluster*fasta) ; do fastANI --ql $fasta --rl $fasta -o ${fasta%.*}.ANI ; done
