#!/usr/bin/env python3.8
from Bio import SeqIO
import sys

multiFasta = sys.argv[1]
outputFile = open(sys.argv[2], "w")

for seq_record in SeqIO.parse(multiFasta, "fasta"):
    outputFile.write(">" + seq_record.id + "\n")
    outputFile.write(str(seq_record.seq.back_transcribe().ungap()) + "\n")

outputFile.close()

# seqRecords = list(SeqIO.parse(multiFasta, "fasta"))
# seq = seqRecords[0]
# print(seq.seq.ungap()[0:50])
# print(seq.seq.back_transcribe().ungap()[0:50])
