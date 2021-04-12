rule species_coverage:
    input:
        f"{dataDir}/{{sample}}/{{sample}}.uniq",
        f"{dataDir}/{{sample}}/dimer{kmerSize}/allOligos_reverse.set.fasta"
    output:
        f"{dataDir}/{{sample}}/evaluation{kmerSize}/allOligos_reverse.set.coverage"
    params:
        splitFilesDir = f"{dataDir}/{{sample}}/splitFiles"
    run:
        import os
        import sys
        import re
        from collections import defaultdict
        from Bio import Entrez
        from Bio import SeqIO
        from Bio.Seq import Seq

        input1 = input[0]
        input2 = input[1]
        output1 = output[0]

        dict_aceIDtaxID=defaultdict(str)
        dict_posiSeq = defaultdict(str)
        dict_posiAce = defaultdict(list)
        dict_taxIDSpecies = defaultdict(str)
        dict_speciesCount = defaultdict(int)

        Entrez.email = "wu.meiju@outlook.com"

        acessionID = os.popen(f"cat {input1}|grep '^>'|cut -f1 -d' '|sed 's/>//g'").read()
        taxID = os.popen(f"cat {input1}|grep '^>'|cut -f1 -d' '|sed 's/>//g'|epost -db nuccore|esummary -db nuccore|xtract -pattern DocumentSummary -element Caption,TaxId|cut -f2").read()
        acessionID = acessionID.split("\n")
        taxID = taxID.split("\n")
        combine = [' '.join(x) for x in zip(acessionID, taxID)]

        for x in combine:
            ls = x.split()
            if len(ls) > 1 :
                ace = x.split()[0]
                tax = x.split()[1]
                dict_aceIDtaxID[ace]=tax

        for aceID in dict_aceIDtaxID:
            taxID = dict_aceIDtaxID[aceID]
            handle = Entrez.efetch(db="taxonomy", id=taxID, mode="text")
            records = Entrez.read(handle)
            for taxon in records:
                for t in taxon["LineageEx"]:
                    if t["Rank"] == "species":
                        species = t["ScientificName"]
                        # print(species)
                        if len(species) > 0 :
                            dict_taxIDSpecies[taxID] = species
                        else:
                            dict_taxIDSpecies[taxID] = "unknow"

        for taxId in dict_taxIDSpecies:
            species = dict_taxIDSpecies[taxId]
            dict_speciesCount[species]+=1

        input2_open = open(input2, "r")
        for line in input2_open.readlines():
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line
                position = header.split("|")[1].split()[1]
                dict_posiSeq[position] = ""
            else:
                kmer = Seq(line)
                oligo = str(kmer.reverse_complement())
                dict_posiSeq[position] = oligo
        input2_open.close()

        for posi in dict_posiSeq:
            seq = dict_posiSeq[posi]
            splitFiles = f"{params.splitFilesDir}/reverse{posi}.fasta"
            accessionID_ls = []
            file = open(splitFiles, "r")
            for line in file:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    accessionID = line.split()[1]
                if not line.startswith(">"):
                    if re.search(seq, line, re.I):
                        accessionID_ls.append(accessionID)
            file.close()
            dict_posiAce[posi]=accessionID_ls

        output1_open = open(output1, "w")
        for posi in dict_posiAce:
            dict=defaultdict(int)
            for acessionID in dict_posiAce[posi]:
                taxID = dict_aceIDtaxID[acessionID]
                species = dict_taxIDSpecies[taxID]
                dict[species]+=1
            for species in dict:
                species_count = dict[species]
                totolCount = dict_speciesCount[species]
                output1_open.write(f"{posi}\t{species}\t{species_count}\t{totolCount}\n")
        output1_open.close()

rule oligoCount:
    input:
        f"{dataDir}/{{sample}}/filtering{kmerSize}/allOligos_reverse.calculated2",
        f"{dataDir}/{{sample}}/filtering{kmerSize}/allOligos_reverse.filtered",
        f"{dataDir}/{{sample}}/checkSpecifity{kmerSize}/allOligos_reverse.filtered.spec.fasta",
        f"{dataDir}/{{sample}}/dimer{kmerSize}/allOligos_reverse.set.fasta"
    output:
        f"{dataDir}/{{sample}}/evaluation{kmerSize}/allOligo_before.tsv",
        f"{dataDir}/{{sample}}/evaluation{kmerSize}/allOligo_after1.tsv",
        f"{dataDir}/{{sample}}/evaluation{kmerSize}/allOligo_after2.tsv",
        f"{dataDir}/{{sample}}/evaluation{kmerSize}/allOligo_after3.tsv"
    run:
        import sys
        import pandas as pd
        from collections import defaultdict

        input1 = input[0]
        input2 = input[1]
        input3 = input[2]
        input4 = input[3]

        tables = {
            input1: output[0],
            input2: output[1]
        }

        fastas = {
            input3: output[2],
            input4: output[3]
        }

        for table in tables:
            df = pd.read_table(table, sep="\t", header=0)
            df = df.sort_values(["position", "kmerCount"], ascending=False)
            df_first = df.groupby("position").first().reset_index()
            df_first=df_first[["position", "kmerCount"]]
            df_first.to_csv(tables[table] , sep="\t", index=False)

        for fasta in fastas:
            dict = defaultdict(list)
            input_open = open(fasta, "r")
            output_open = open(fastas[fasta], "w")

            for line in input_open:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    ls = line.split("|")
                    position = ls[1].split()[1]
                    count = int(ls[2].split()[1])
                    # print(ls)
                    # print(position, count)
                    dict[position].append(count)

            output_open.write(f"position\tkmerCount\n")
            for key in dict:
                ls = sorted(dict[key])
                output_open.write(f"{key}\t{ls[-1]}\n")

            input_open.close()
            output_open.close()
