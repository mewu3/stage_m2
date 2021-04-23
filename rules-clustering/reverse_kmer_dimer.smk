rule allHeterodimerCheck:
    input:
        f"{dataDir}/{{sample}}{clustering_identity}/checkSpecifity{kmerSize}/allOligos_reverse.filtered.spec.fasta"
    output:
        f"{dataDir}/{{sample}}{clustering_identity}/dimer{kmerSize}/allOligos_reverse.filtered.spec.heterodimer"
    params:
        monovalentConc = config["oligotm"]["monovalent-conc"],
        divalentConc = config["oligotm"]["divalent-conc"],
        dNTPConc = config["oligotm"]["dNTP-conc"],
        dnaConc = config["oligotm"]["dna-conc"],
        thermodynamicPara = config["oligotm"]["thermodynamic-para"],
        saltCorrelation = config["oligotm"]["salt-correlation"]
    run:
        import os
        import primer3
        from collections import defaultdict

        inputOpen = open(input[0], "r")
        outputOpen = open(output[0], "w")

        dict_positionIdKmer = defaultdict(lambda: defaultdict(str))

        for line in inputOpen.readlines():
            line = line.rstrip("\n")
            if line.startswith(">"):
                ID = line.split("|")[0].lstrip(">p").strip()
                position = line.split("|")[1].split()[1]
                dict_positionIdKmer[position][ID]=""
            else:
                dict_positionIdKmer[position][ID]=line

        outputOpen.write(f"position1\tID1\tposition2\tID2\theterodimer\n")
        for position1 in dict_positionIdKmer:
            for position2 in dict_positionIdKmer:
                if position1 != position2:
                    for ID1 in dict_positionIdKmer[position1]:
                        for ID2 in dict_positionIdKmer[position2]:
                            k1=dict_positionIdKmer[position1][ID1]
                            k2=dict_positionIdKmer[position2][ID2]
                            heterodimer = primer3.calcHeterodimer(k1, k2,
                                                                  mv_conc = params.monovalentConc,
                                                                  dv_conc = params.divalentConc,
                                                                  dntp_conc = params.dNTPConc,
                                                                  dna_conc = params.dnaConc).dg
                            if float(heterodimer) > -9000:
                                outputOpen.write(f"{position1}\t{ID1}\t{position2}\t{ID2}\t{heterodimer}\n")
                            # else:
                            #     list.append("1")

        # formatedHeader = "\t".join(map(str,ID_list))
        # outputOpen.write(f"ID\t{formatedHeader}\n")

        # for id1 in dictIDKmer:
        #     k1 = dictIDKmer[id1]
        #     list = []
        #     for id2 in dictIDKmer:
        #         k2 = dictIDKmer[id2]
        #         heterodimer = primer3.calcHeterodimer(k1, k2,
        #                                               mv_conc = params.monovalentConc,
        #                                               dv_conc = params.divalentConc,
        #                                               dntp_conc = params.dNTPConc,
        #                                               dna_conc = params.dnaConc).dg
        #         if float(heterodimer) > -9000:
        #             list.append("0")
        #         else:
        #             list.append("1")
        #     formatedList = "\t".join(map(str,list))
            # outputOpen.write(f"{id1}\t{formatedList}\n")

        inputOpen.close()
        outputOpen.close()

rule oligoSet:
    input:
        f"{dataDir}/{{sample}}/dimer{kmerSize}/allOligos_reverse.filtered.spec.heterodimer",
        f"{dataDir}/{{sample}}/checkSpecifity{kmerSize}/allOligos_reverse.filtered.spec.fasta"
    output:
        f"{dataDir}/{{sample}}/dimer{kmerSize}/allOligos_reverse.set",
        f"{dataDir}/{{sample}}/dimer{kmerSize}/allOligos_reverse.set.fasta"
    params:
        monovalentConc = config["oligotm"]["monovalent-conc"],
        divalentConc = config["oligotm"]["divalent-conc"],
        dNTPConc = config["oligotm"]["dNTP-conc"],
        dnaConc = config["oligotm"]["dna-conc"],
        thermodynamicPara = config["oligotm"]["thermodynamic-para"],
        saltCorrelation = config["oligotm"]["salt-correlation"]
    run:
        from collections import defaultdict
        import re
        import pandas as pd
        import primer3
        from Bio import SeqIO

        input1_open = open(input[0], "r")
        # input2_open = open(input[1], "r")
        output1_open = open(output[0], "w")
        output2_open = open(output[1], "w")

        dict = defaultdict(lambda: defaultdict(lambda:defaultdict(list)))



        for line in input1_open.readlines()[1:]:
            line = line.rstrip("\n")
            ls = line.split('\t')
            posi1 = str(ls[0])
            oligo1=ls[1]
            posi2=str(ls[2])
            oligo2=ls[3]
            if posi2 not in dict and oligo2 not in dict[posi1]:
                dict[posi1][oligo1][posi2].append(oligo2)

        oligoSet = set()

        keyCount = len(list(dict.keys()))

        for n in range(0, keyCount-1):

            firstPosi1 = list(dict.keys())[n]
            firstOligo1 = list(dict[firstPosi1].keys())[0]
            firstPosi2 = list(dict[firstPosi1][firstOligo1].keys())[0]
            firstOligo2 = dict[firstPosi1][firstOligo1][firstPosi2][0]
            # print(firstPosi1, firstOligo1, firstPosi2, firstOligo2)
            oligoSet.add(firstOligo1)

            secondPosi1 = list(dict.keys())[n+1]
            secondOligo1 = firstOligo2
            secondPosi2 = list(dict[secondPosi1][secondOligo1].keys())[0]
            secondOligo2 = dict[secondPosi1][secondOligo1][secondPosi2][0]
            # print(secondPosi1, secondOligo1, secondPosi2, secondOligo2)
            oligoSet.add(secondOligo2)

        record_dict = SeqIO.to_dict(SeqIO.parse(input[1], "fasta"))

        oligoSet = list(map(lambda x: "p"+x, oligoSet))
        formatedHeader = "\t".join(map(str,oligoSet))
        output1_open.write(f"ID\t{formatedHeader}\n")

        for id1 in oligoSet:
            ls=[]
            seq1 = str(record_dict[id1].seq)
            output2_open.write(record_dict[id1].format("fasta"))
            for id2 in oligoSet:
                seq2 = str(record_dict[id2].seq)
                heterodimer = primer3.calcHeterodimer(seq1, seq2,
                                                      mv_conc = params.monovalentConc,
                                                      dv_conc = params.divalentConc,
                                                      dntp_conc = params.dNTPConc,
                                                      dna_conc = params.dnaConc).dg
                if float(heterodimer) > -9000:
                    ls.append("0")
                else:
                    ls.append("1")
                    print(f"{id1} is not compatible with {id2}")
            formatedList = "\t".join(map(str,ls))
            output1_open.write(f"{id1}\t{formatedList}\n")


        # formatedHeader = "\t".join(map(str,oligoSet))
        # output1_open.write(f"ID\t{formatedHeader}\n")
        # for o1 in oligoSet:
        #     ls=[]
        #     for o2 in oligoSet:
        #         heterodimer = primer3.calcHeterodimer(o1, o2,
        #                                               mv_conc = params.monovalentConc,
        #                                               dv_conc = params.divalentConc,
        #                                               dntp_conc = params.dNTPConc,
        #                                               dna_conc = params.dnaConc).dg
        #         if float(heterodimer) > -9000:
        #             ls.append("0")
        #         else:
        #             ls.append("1")
        #             print(f"{o1} is not compatible with {o2}")
        #     formatedList = "\t".join(map(str,ls))
        #     output1_open.write(f"{o1}\t{formatedList}\n")

        input1_open.close()
        # input2_open.close()
        output1_open.close()
        output2_open.close()
