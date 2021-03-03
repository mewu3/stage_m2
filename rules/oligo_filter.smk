# kmerCounter = "kmc"

rule calculte_tm:
    input:
        datadir + "/kmerCounting/" + kmerCounter + "/{prefix}.sort"
    output:
        datadir + "/filter/" + kmerCounter + "/{prefix}.txt"
    params:
        monovalentConc = config["oligotm"]["monovalent-conc"],
        divalentConc = config["oligotm"]["divalent-conc"],
        dNTPConc = config["oligotm"]["dNTP-conc"],
        dnaConc = config["oligotm"]["dna-conc"],
        thermodynamicPara = config["oligotm"]["thermodynamic-para"],
        saltCorrelation = config["oligotm"]["salt-correlation"]
    run:
        import os
        import subprocess

        output = open(output[0], "w")
        output.write("{} {} {}\n".format("kmer", "kmer_count", "Tm"))

        with open(input[0], "r") as input:
            for li in input:
                li = li.rstrip("\n")
                ls = li.split()
                cmd = "lib/primer3/src/oligotm -mv {} -dv {} -n {} -d {} -tp {} -sc {} {}".format(params.monovalentConc, params.divalentConc, params.dNTPConc, params.dnaConc, params.thermodynamicPara, params.saltCorrelation, ls[0])
                tm = os.popen(cmd).read()
                output.write("{} {} {}".format(ls[0], ls[1], tm))

        output.close()

rule filterOut_Tm:
    input:
        datadir + "/filter/" + kmerCounter + "/{prefix}.txt"
    output:
        datadir + "/filter/" + kmerCounter + "/{prefix}.table"
    run:
        import pandas as pd

        df = pd.read_table(input[0], sep=" ", header=0)

        df_filtered = df[(df["Tm"] >= 60) & (df["Tm"] <= 70)]

        df_filtered.to_csv(output[0], sep='\t')
