rule forward_TmCalculate:
    input:
        datadir + "/{sample}/dsk/forward_{seg}.kCountSorted",
    output:
        datadir + "/{sample}/primer3_Tm/forward_{seg}.Tm",
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
        output.write("{}\t{}\t{}\t{}\n".format("kmer", "kmer_count", "position", "Tm"))

        with open(input[0], "r") as input:
            for li in input:
                li = li.rstrip("\n")
                ls = li.split()
                cmd = "lib/primer3/src/oligotm -mv {} -dv {} -n {} -d {} -tp {} -sc {} {}".format(params.monovalentConc, params.divalentConc, params.dNTPConc, params.dnaConc, params.thermodynamicPara, params.saltCorrelation, ls[0])
                tm = os.popen(cmd).read()
                output.write("{}\t{}\t{}\t{}".format(ls[0], ls[1], wildcards.seg, tm))

        output.close()

checkpoint filterOut_forward:
    input:
        datadir + "/{sample}/primer3_Tm/forward_{seg}.Tm"
    output:
        datadir + "/{sample}/primer3_Tm/forward_{seg}.TmFilterered"
    run:
        import pandas as pd
        import os

        df = pd.read_table(input[0], sep="\t", header=0)
        df_filtered = df[(df["Tm"] >= 60) & (df["Tm"] <= 70)]
        df_filtered.to_csv(output[0], sep='\t', index=False)

def aggregate_input(wildcards):
    checkpoint_output = datadir + "/{sample}/primer3_Tm/"
    return expand(datadir + "/{sample}/primer3_Tm/forward_{seg}.TmFilterered",
                  sample = wildcards.sample,
                  seg = glob_wildcards(os.path.join(checkpoint_output, "forward_{seg}.TmFilterered")).seg)

rule merge_all_kmers:
    input:
        aggregate_input
    output:
        datadir + "/{sample}/primer3_Tm/all_oligo"
    shell:
        "sed -s 1d {input} >> {output}"

rule transforme_to_fasta:
    input:
        datadir + "/{sample}/primer3_Tm/all_oligo"
    output:
        datadir + "/{sample}/primer3_Tm/all_oligo.fasta"
    run:
        import os

        outFile = open(output[0], "w")

        with open(input[0], "r") as file:
            for li in file:
                li = li.rstrip("\n")
                ls = li.split()
                kmer = ls[0]
                count = ls[1]
                position = ls[2]
                tm = ls[3]
                outFile.write(">position {} |count {} |Tm {}\n {}\n".format(position, count, tm, kmer))

        outFile.close()

rule index_fasta:
    input:
        datadir + "/{sample}/primer3_Tm/all_oligo.fasta"
    output:
        datadir + "/{sample}/primer3_Tm/all_oligo.fasta.fai"
        datadir + "/{sample}/primer3_Tm/all_oligo.fasta.json"
        datadir + "/{sample}/primer3_Tm/all_oligo.fasta.log"
        datadir + "/{sample}/primer3_Tm/all_oligo.fasta.primerqc"
        datadir + "/{sample}/primer3_Tm/all_oligo.fasta.primerqc.fai"
    shell:
        "lib/mfeprimer-3.2.1-linux-amd64 -i {input}"
