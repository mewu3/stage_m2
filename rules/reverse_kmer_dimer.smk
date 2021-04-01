rule allHeterodimerCheck:
    input:
        f"{dataDir}/{{sample}}/checkSpecifity/allOligos_reverse.filtered.spec.fasta"
    output:
        f"{dataDir}/{{sample}}/checkSpecifity/allOligos_reverse.filtered.spec.heterodimer"
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
        f"{dataDir}/{{sample}}/checkSpecifity/allOligos_reverse.filtered.spec.heterodimer"
    output:
        f"{dataDir}/{{sample}}/checkSpecifity/allOligos_reverse.set"
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

        inputOpen = open(input[0], "r")
        outputOpen = open(output[0], "w")

        dict = defaultdict(lambda: defaultdict(lambda:defaultdict(list)))
        # dict=defaultdict(list)
        # seenPosition1 = []
        # seenPosition2 =[]
        # seenOligo = []
        # seenOligo1 = []
        # seenOligo2 = []


        for line in inputOpen.readlines()[1:]:
            line = line.rstrip("\n")
            ls = line.split('\t')
            posi1 = str(ls[0])
            oligo1=ls[1]
            posi2=str(ls[2])
            oligo2=ls[3]
            # if posi1 not in dict:
            #     dict[posi1]={}
            # if posi2 not in dict and oligo2 not in dict[posi1][posi2][oligo1]:
            #     dict[posi1][posi2]={}
            if posi2 not in dict and oligo2 not in dict[posi1]:
                dict[posi1][oligo1][posi2].append(oligo2)
            # if posi1 not in dict:
            #     dict[posi1]={}
                # seenPosition1.append(o)
            # if posi2 not in dict[posi1] and posi2 not in dict:
                # dict[posi1][posi2]=[oligo1, oligo2]
                # dict[posi1][posi2]={}
                # dict[posi1][posi2][oligo1]=oligo2
        #     if posi2 not in seenPosition2 :
        #         seenPosition1.append(posi2)
        #         print(line)
                # seenPosition.append(posi2)
                # if oligo1 not in seenOligo and oligo2 not in seenOligo:
                #     seenOligo.append(oligo1)
                #     seenOligo.append(oligo2)
                #     print(line)
                # if oligo1 in seenOligo or oligo2 in seenOligo:
                #     oligoSet.add(oligo1)
                #     oligoSet.add(oligo2)
                #     print(line)
                # if posi2 not in dict[posi1] and posi2 != posi1 :
                #         dict[posi1][posi2] = [oligo1, oligo2]
                #         print(posi1, posi2, oligo1, oligo2)

        # print(dict)
        # print(oligoSet)

        oligoSet = []

        # firstPosi1 = list(dict.keys())[0]
        # firstOligo1 = list(dict[firstPosi1].keys())[0]
        # firstPosi2 = list(dict[firstPosi1][firstOligo1].keys())[0]
        # firstOligo2 = dict[firstPosi1][firstOligo1][firstPosi2][0]
        # print(firstOligo1, firstOligo2)

        # secondPosi1 = list(dict.keys())[1]
        # secondOligo1 = firstOligo2
        # secondPosi2 = list(dict[secondPosi1][secondOligo1].keys())[0]
        # secondOligo2 = dict[secondPosi1][secondOligo1][secondPosi2][0]
        # print(secondOligo1, secondOligo2)

        keyCount = len(list(dict.keys()))

        for n in range(0, keyCount-1):

            firstPosi1 = list(dict.keys())[n]
            firstOligo1 = list(dict[firstPosi1].keys())[0]
            firstPosi2 = list(dict[firstPosi1][firstOligo1].keys())[0]
            firstOligo2 = dict[firstPosi1][firstOligo1][firstPosi2][0]
            # print(firstPosi1, firstOligo1, firstPosi2, firstOligo2)
            oligoSet.append(firstOligo1)

            secondPosi1 = list(dict.keys())[n+1]
            secondOligo1 = firstOligo2
            secondPosi2 = list(dict[secondPosi1][secondOligo1].keys())[0]
            secondOligo2 = dict[secondPosi1][secondOligo1][secondPosi2][0]
            # print(secondPosi1, secondOligo1, secondPosi2, secondOligo2)
            oligoSet.append(secondOligo2)

            # if secondOligo2 in dict[firstPosi1][firstOligo1][secondPosi2]:
            #     print(firstOligo1, secondOligo2)

        formatedHeader = "\t".join(map(str,oligoSet))
        outputOpen.write(f"ID\t{formatedHeader}\n")
        for o1 in oligoSet:
            ls=[]
            for o2 in oligoSet:
                heterodimer = primer3.calcHeterodimer(o1, o2,
                                                      mv_conc = params.monovalentConc,
                                                      dv_conc = params.divalentConc,
                                                      dntp_conc = params.dNTPConc,
                                                      dna_conc = params.dnaConc).dg
                if float(heterodimer) > -9000:
                    ls.append("0")
                else:
                    ls.append("1")
                    print(f"{o1} is not compatible with {o2}")
            formatedList = "\t".join(map(str,ls))
            outputOpen.write(f"{o1}\t{formatedList}\n")

        # oligoSet.append(firstOligo1)
        # oligoSet.append(firstOligo2)
        #
        # for p1 in dict:
        #     for k1 in dict[p1]:
        #         for p2 in dict[p1][k1]:
        #             if dict[p1][k1][p2]



                # oligo1 = dict[p1][p2][0]
                # oligo2 = dict[p1][p2][1]
                # print(p1, p2, oligo1, oligo2)
                # print(dict[p1][p2])
                # if oligo1 not in seenOligo1:
                #     seenOligo1.append(oligo1)
                # if oligo2 not in seenOligo2:
                #     seenOligo2.append(oligo2)
                # if oligo1 in seenOligo2:
                #     print(dict[p1][p2])
                # print(seenOligo1)
                # print(seenOligo2)
                # if oligo1 in seenOligo2:
                #     print(dict[p1][p2])
                #     if oligo1 in seenOligo or oligo2 in seenOligo:
                #         oligoSet.append(dict[p1][p2])
                #         print(dict[p1][p2])


        # df = pd.read_table(input[0], sep="\t", header=0)
        # df_group = df.groupby(["position1","position2"])
        # print(df.groupby("position1")["position1"].count())
        #
        # pd.set_option("display.max_rows", None, "display.max_columns", None)

        # df_first = df_group.first()
        # print(df_first["ID1"].unique())
        # print(df_first["ID2"].unique())
        # print(df_first)
        # ~ print(df_group.apply(list).to_dict())

        # ~ for k, v in df_group:
            # ~ print(v.first())

        inputOpen.close()
        outputOpen.close()
