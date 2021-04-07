import os
from collections import defaultdict
from Bio import Entrez

uniq_seq = "/datater/wu/data/enterovirus/enterovirus.uniq"
uniq_seq_DB = "/datater/wu/data/enterovirus/enterovirus.uniq.DB"
oligoSet = "/datater/wu/data/enterovirus/checkSpecifity/allOligos_reverse.set.fasta"
oligoSet_blastn = "/datater/wu/data/enterovirus/checkSpecifity/allOligos_reverse.set.fasta.coverage"
taxID_in="/home/meijun/Téléchargements/taxid_input.txt"
taxID_out = "/home/meijun/Téléchargements/taxid_out.txt"
coverage_tax = "/home/meijun/Téléchargements/taxCoverage.txt"
execute = False

if execute:
    os.system(f"makeblastdb -in {uniq_seq} -out {uniq_seq_DB} -dbtype nucl -parse_seqids")
    os.system(f"blastn -task blastn-short -db {uniq_seq_DB} -word_size 7 -max_target_seqs 5 -query {oligoSet} -out {oligoSet_blastn} -outfmt '6 qacc sacc stitle pident length mismatch gapopen sstart send evalue'")
    os.system(f"cat {oligoSet_blastn} |cut -f2 |epost -db nuccore |esummary -db nuccore |xtract -pattern DocumentSummary -element Caption,TaxId |cut -f2")
    os.system(f"egrep '^>' {uniq_seq} |cut -f1 -d" " |sed 's/>//g' |epost -db nuccore |esummary -db nuccore |xtract -pattern DocumentSummary -element Caption,TaxId |cut -f2")

taxID_in_open = open(taxID_in, "r")
taxID_out_open = open(taxID_out, "r")
coverage_tax_open = open(coverage_tax, "w")

def count_taxid(file):
    dict_count = defaultdict(int)
    for line in file:
        line = line.rstrip("\n")
        dict_count[line] += 1
    return dict_count

dict_in = count_taxid(taxID_in_open)
dict_out = count_taxid(taxID_out_open)

Entrez.email = "wu.meiju@outlook.com"

out_count_sum = 0
in_count_sum = 0

dict_specie_inCount = defaultdict(int)
dict_specie_outCount = defaultdict(int)

for id in dict_in:
    handle = Entrez.efetch(db="taxonomy", id=id, mode="text")
    records = Entrez.read(handle)
    species = "unknow"
    for taxon in records:
        for t in taxon["LineageEx"]:
            if t["Rank"] == "species":
                species = t["ScientificName"]
    # dict_specie_inCount[species] = 0
    # dict_specie_outCount[species] = 0
    in_count_sum += dict_in[id]
    ### variants levels
    if id in dict_out:
        percentage = dict_out[id]/dict_in[id]
        out_count = dict_out[id]
        in_count = dict_in[id]
        out_count_sum += dict_out[id]
        dict_specie_inCount[species] += in_count
        dict_specie_outCount[species] += out_count
        coverage_tax_open.write(f"{species}\t{id}\t{out_count}\t{in_count}\t{percentage}\n")
    else:
        percentage = 0
        out_count = 0
        in_count = dict_in[id]
        dict_specie_inCount[species] += in_count
        dict_specie_outCount[species] += out_count
        coverage_tax_open.write(f"{species}\t{id}\t{out_count}\t{in_count}\t{percentage}\n")

CoveragePercentage = out_count_sum/in_count_sum *100
coverage_tax_open.write(f"The output oligoSet target (variant levels) {out_count_sum} out of {in_count_sum}: {CoveragePercentage:.3f}%.\n")

for species in dict_specie_inCount:
    inCount = dict_specie_inCount[species]
    out_count = dict_specie_outCount[species]
    percentage = out_count/inCount
    coverage_tax_open.write(f"The output oligoSet target {species} : {out_count} / {inCount} = {percentage:.3f} \n")

taxID_in_open.close()
taxID_out_open.close()
coverage_tax_open.close()
