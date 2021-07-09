#!/usr/bin/env python3.8

import sys, os
from Bio import SeqIO
from pathlib import Path
import pandas as pd

def splitMultiFasta(multiFasta, outputDir):
    Path(outputDir).mkdir(parents=True, exist_ok=True)
    with open(multiFasta) as multiFasta:
        records = SeqIO.parse(multiFasta, "fasta")
        for record in records:
            outputFile = f"{outputDir}/{record.id}.fasta"
            with open(outputFile, "w") as outputFile:
                SeqIO.write(record, outputFile, "fasta")
