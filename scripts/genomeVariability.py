from Bio import SeqIO
import sys
import pandas as pd
import os
import numpy as np
import fnmatch
import glob
import shutil

## split multifasta to individual fasta ### 
# f_open = open(sys.argv[1], "r")
# for rec in SeqIO.parse(f_open, "fasta"):
#    id = rec.id
#    seq = rec.seq
#    id_file = open(id, "w")
#    id_file.write(">"+str(id)+"\n"+str(seq))
#    id_file.close()
# f_open.close()

### move individual fasta to specie folder ### 
# df = pd.read_table(
#     "/datater/wu/MSSPE/scripts/test.tsv",
#     sep = "\t",
#     header = 0,
#     index_col = 0
# )
#
# dict_taxID_aceID = {k: g["accession"].tolist() for k, g in df[["tax_id", "accession"]].groupby("tax_id")}
# dict_spe_taxID = {k: g["tax_id"].tolist() for k, g in df[["tax_id", "species"]].groupby("species")}
#
# dict_spe_taxID_aceID = {}
#
# for spe in dict_spe_taxID:
#     dict_spe_taxID_aceID[spe]={}
#     for taxID in dict_spe_taxID[spe]:
#         dict_spe_taxID_aceID[spe][taxID] = dict_taxID_aceID[taxID]

# for spe in dict_spe_taxID_aceID:
#     for taxID in dict_spe_taxID_aceID[spe]:
#         speName = spe.replace(" ", "_")
#         path = f"{speName}/{taxID}"
#         # os.makedirs(path, exist_ok=True)
#         for aceID in dict_spe_taxID_aceID[spe][taxID]:
#             commande = f"mv {aceID}* {path}"
#             os.system(commande)

### create fastANI input files for each taxID ### 
# ls_species = glob.glob(f"/datater/wu/data/coronavirus/*/")
#
# for specie in ls_species:
#     for taxID in glob.glob(f"{specie}*[!txt]"):
#         taxid = taxID.split("/")[-1]
#         fastANI_input = f"{taxID}/{taxid}.txt"
#         # print(fastANI_input)
#         fastANI_input_open = open(fastANI_input, "w")
#         files = glob.glob(f"{taxID}/*[!txt]")
#         # print(files)
#         for f in files:
#             fastANI_input_open.write(f"{f}\n")
#         fastANI_input_open.close()

### crate fastANI input files for each species ### 
# ls_species = glob.glob(f"/datater/wu/data/coronavirus/*/")
#
# for specie in ls_species:
#     specieName = specie.split("/")[-2]
#     specieFile = f"{specie}{specieName}.txt"
#     specieFile = open(specieFile, "w")
#     for taxID in glob.glob(f"{specie}*[!txt]"):
#         taxid = taxID.split("/")[-1]
#         fastANI_input = f"{taxID}/{taxid}.txt"
#         fastANI_input = open(fastANI_input, "r")
#         for line in fastANI_input:
#             specieFile.write(line)
#         fastANI_input.close()
#     specieFile.close()

### all and all ### 
# ls_species = glob.glob(f"/datater/wu/data/enterovirus/*/")
# all_files = f"/datater/wu/data/enterovirus/enterovirus.txt"
# all_files = open(all_files, "w")
# for specie in ls_species:
#     specieName = specie.split("/")[-2]
#     specieFile = f"{specie}{specieName}.txt"
#     specieFile = open(specieFile, "r")
#     for li in specieFile:
#         all_files.write(li)
#     specieFile.close()
# all_files.close()
