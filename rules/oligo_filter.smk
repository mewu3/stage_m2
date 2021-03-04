rule dskTm_forward_calculate:
    input:
        datadir + "/kmerCounting/{sample}/forward/{seg}_dsk.sort",
    output:
        datadir + "/filtering/{sample}/forward/{seg}_dsk.table",
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
                output.write("{} {} {}\n".format(ls[0], ls[1], tm))

        output.close()

rule filterOut_forward_dskTm:
    input:
        datadir + "/filtering/{sample}/forward/{seg}_dsk.table"
    output:
        datadir + "/filtering/{sample}/forward/{seg}_dsk.tsv"
    run:
        import pandas as pd
        import os

        df = pd.read_table(input[0], sep=" ", header=0)
        df_filtered = df[(df["Tm"] >= 60) & (df["Tm"] <= 70)]
        df_filtered.to_csv(output[0], sep='\t')


rule dskTm_reverse_calculate:
    input:
        datadir + "/kmerCounting/{sample}/reverse/{seg}_dsk.sort",
    output:
        datadir + "/filtering/{sample}/reverse/{seg}_dsk.table",
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
                output.write("{} {} {}\n".format(ls[0], ls[1], tm))

        output.close()

rule filterOut_reverse_dskTm:
    input:
        datadir + "/filtering/{sample}/reverse/{seg}_dsk.table"
    output:
        datadir + "/filtering/{sample}/reverse/{seg}_dsk.tsv"
    run:
        import pandas as pd
        import os

        df = pd.read_table(input[0], sep=" ", header=0)
        df_filtered = df[(df["Tm"] >= 60) & (df["Tm"] <= 70)]
        df_filtered.to_csv(output[0], sep='\t')
