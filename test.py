import glob
import re

path = "/datater/wu/data/splitFiles"
samples = ["enterovirus"]

for file in glob.glob(path+"*.fasta", recursive=True):
    print(file)
