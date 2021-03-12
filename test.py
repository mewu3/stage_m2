# input1 = open("/datater/wu/data/enterovirus/filtering/forward_2032-2536.nessieOut", "r")
input2 = open("/datater/wu/data/enterovirus/filtering/forward_2032-2536.fasta", "r")
outputFile = open("/datater/wu/data/enterovirus/filtering/forward_2032-2536.calculated2", "w")

dict_id_lc = {}

with open("/datater/wu/data/enterovirus/filtering/forward_2032-2536.nessieOut", "r") as input1 :
    next(input1)
    for li1 in input1:
        li1 = li1.rstrip("\n")
        if li1.startswith(">"):
            id = li1.split("|")[0].lstrip(">")
            dict_id_lc[id] = 0
        else:
            LC = float(li1.split(":")[1])
            dict_id_lc[id] = LC

for li2 in input2:
    li2 = li2.rstrip("\n")
    # id_ = ""
    # position = ""
    # count = ""
    # cg = ""
    # tm = ""
    # homodimer_dG = ""
    # hairpin_dG = ""
    # lc = ""
    # seq = ""
    list_ = []

    if li2.startswith(">"):
        ls = li2.split("|")
        id_ = ls[0].lstrip(">")
        position = ls[1].split()[1]
        count = ls[2].split()[1]
        cg = ls[3].split()[1]
        tm = ls[4].split()[1]
        homodimer_dG = ls[5].split()[1]
        hairpin_dG = ls[6].split()[1]
        outputFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(id_, position, count, cg, tm, homodimer_dG, hairpin_dG))
        for id in dict_id_lc:
            if id_ == id :
                lc = str(dict_id_lc[id_])
                outputFile.write("{}\t".format(lc))
    else:
        seq = li2
        outputFile.write("{}\n".format(seq))
    # print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(id_, seq, position, count, cg, tm, homodimer_dG, hairpin_dG))
