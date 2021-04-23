rule get taxID_and_specieName:
    input:
        f"{dataDir}/{{sample}}/{{sample}}.uniq"
    output:
        f"{dataDir}/{{sample}}/evaluation{kmerSize}/aceID-taxID-species.tsv"
    run:
        import os
        import sys
        import re
        from collections import defaultdict
        from Bio import Entrez
        from Bio import SeqIO
        from Bio.Seq import Seq

        input1 = input[0]
        output1 = output[0]

        dict_aceIDtaxID=defaultdict(str)

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

        # handle = Entrez.efetch(db="taxonomy", id="2760819", mode="text")
        # records = Entrez.read(handle)
        # for taxon in records:
        #     print(taxon)
            # if taxon["Rank"] == "species":
            #     species = taxon["ScientificName"]

        # handle = Entrez.efetch(db="taxonomy", id="86107", mode="text")
        # records = Entrez.read(handle)
        # for taxon in records:
        #     if taxon["Rank"] == "serotype":
        #         for t in taxon["LineageEx"]:
        #             if t["Rank"] == "species":
        #                 species = t["ScientificName"]
                        # output1_open.write(f"{aceID}\t{taxID}\t{species}\n")

        output1_open = open(output1, "w")
        for aceID in dict_aceIDtaxID:
            taxID = dict_aceIDtaxID[aceID]
            handle = Entrez.efetch(db="taxonomy", id=taxID, mode="text")
            records = Entrez.read(handle)
            species = "unknow"
            for taxon in records:
                if taxon["Rank"] == "species":
                    species = taxon["ScientificName"]
                # if taxon["Rank"] == "serotype":
                #     for t in taxon["LineageEx"]:
                #         if t["Rank"] == "species":
                #             species = t["ScientificName"]
                # if taxon["Rank"] == "no rank":
                else:
                    for t in taxon["LineageEx"]:
                        if t["Rank"] == "species":
                            species = t["ScientificName"]
                output1_open.write(f"{aceID}\t{taxID}\t{species}\n")
        output1_open.close()

rule species_coverage:
    input:
        f"{dataDir}/{{sample}}/dimer{kmerSize}/allOligos_reverse.set.fasta",
        f"{dataDir}/{{sample}}/evaluation{kmerSize}/aceID-taxID-species.tsv"
    output:
        f"{dataDir}/{{sample}}/evaluation{kmerSize}/allOligos_reverse.set.coverage"
    params:
        splitFilesDir = f"{dataDir}/{{sample}}/splitFiles"
    run:
        import os
        import sys
        import re
        from collections import defaultdict
        from Bio import SeqIO
        from Bio.Seq import Seq

        input1 = input[0]
        input2 = input[1]
        output1 = output[0]

        dict_posiSeq = defaultdict(str)
        dict_posiAce = defaultdict(list)
        dict_aceIdSpecie = defaultdict(str)
        dict_speciesCount = defaultdict(int)

        input1_open = open(input1, "r")
        for line in input1_open.readlines():
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line
                position = header.split("|")[1].split()[1]
                dict_posiSeq[position] = ""
            else:
                kmer = Seq(line)
                oligo = str(kmer.reverse_complement())
                dict_posiSeq[position] = oligo
        input1_open.close()

        for posi in dict_posiSeq:
            seq = dict_posiSeq[posi]
            splitFiles = f"{params.splitFilesDir}/reverse{posi}.fasta"
            accessionID_ls = []
            record_dict = SeqIO.to_dict(SeqIO.parse(splitFiles, "fasta"))
            for record in record_dict:
                if re.search(seq, str(record_dict[record].seq), re.I):
                    accessionID_ls.append(record_dict[record].id)
            dict_posiAce[posi]=accessionID_ls

        input2_open = open(input2, "r")
        for li in input2_open:
            li = li.rstrip("\n")
            ls = li.split("\t")
            aceId = ls[0]
            specie = ls[2]
            dict_aceIdSpecie[aceId]=specie
            dict_speciesCount[specie] += 1
        input2_open.close()

        output1_open = open(output1, "w")
        for posi in dict_posiAce:
            dict=defaultdict(int)
            for acessionID in dict_posiAce[posi]:
                spec = dict_aceIdSpecie[acessionID]
                dict[spec]+=1
            for specie in dict:
                species_count = dict[specie]
                totolCount = dict_speciesCount[specie]
                output1_open.write(f"{posi}\t{specie}\t{species_count}\t{totolCount}\n")
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
