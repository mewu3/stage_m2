rule forward_TmCalculate:
    input:
        datadir + "/{sample}/" + kmerCounter + "/forward_{seg}.kCountSorted"
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

# rule transforme_to_fasta:
#     input:
#         datadir + "/{sample}/primer3_Tm/all_oligo"
#     output:
#         datadir + "/{sample}/primer3_Tm/all_oligo.fasta"
#     run:
#         import os
#
#         outFile = open(output[0], "w")
#
#         n = 0
#         with open(input[0], "r") as file:
#             for li in file:
#                 li = li.rstrip("\n")
#                 ls = li.split()
#                 kmer = ls[0]
#                 count = ls[1]
#                 position = ls[2]
#                 tm = ls[3]
#                 outFile.write(">p{} |position {} |count {} |Tm {}\n {}\n".format(n, position, count, tm, kmer))
#                 n += 1
#
#         outFile.close()

# rule index_fasta:
#     input:
#         datadir + "/{sample}/primer3_Tm/all_oligo.fasta"
#     output:
#         expand(datadir + "/{{sample}}/primer3_Tm/all_oligo.fasta.{ext}",
#             ext = ["fai", "json", "log", "primerqc", "primerqc.fai"]
#         )
#     shell:
#         "lib/mfeprimer-3.2.1-linux-amd64 index -i {input}"

# rule hairpin_mfeprimer:
#     input:
#         datadir + "/{sample}/primer3_Tm/all_oligo.fasta"
#     output:
#         datadir + "/{sample}/hairpin_out"
#     params:
#         deltaG = config["mfeprimer"]["deltaG-cutoff"],
#         Tm = config["mfeprimer"]["Tm-cutoff"],
#         monovalent = config["mfeprimer"]["monovalent-conc"],
#         divalent = config["mfeprimer"]["divalent-conc"],
#         NPT = config["mfeprimer"]["dNTP-conc"],
#         oligo = config["mfeprimer"]["oligo-conc"]
#     shell:
#         "lib/mfeprimer-3.2.1-linux-amd64 hairpin -i {input} -o {output} \
#         -t {params.Tm} \
#         -d {params.deltaG} \
#         --mono {params.monovalent} \
#         --diva {params.divalent} \
#         --dntp {params.NPT} \
#         --oligo {params.oligo}"

rule hairpin_ntthal:
    input:
        datadir + "/{sample}/primer3_Tm/all_oligo"
    output:
        datadir + "/{sample}/primer3_ntthal/allOligo_hairpin"
    run:
        outFile = open(output[0], "w")

        with open(input[0], "r") as file:
            for li in file:
                li = li.rstrip("\n")
                ls = li.split()
                kmer = ls[0]
                stdout = os.popen("lib/primer3/src/ntthal -a HAIRPIN -s1 {}".format(kmer)).read()
                list = stdout.split()
                if len(list) > 6 :
                    deltaG = float(list[14])
                    if deltaG > -9000 :
                        outFile.write("{} {}\n".format(li, deltaG))
        outFile.close()

rule dimer_ntthal:
    input:
        datadir + "/{sample}/primer3_ntthal/allOligo_hairpin"
    output:
        datadir + "{sample}/ptrimer3_ntthal/allOligo_hairpin_dimer"
    run:
        kmerList = []

        # outFile = open(output[0], "w")

        with open(input[0], "r") as file:
            for li in file:
                    li = li.rstrip("\n")
                    ls = li.split()
                    kmerList.append(ls[0])

        for k1 in kmerList:
            for k2 in kmerList:
                if k1 != k2:
                    stdout = os.popen("lib/primer3/src/ntthal -s1 {} -s2 {}".format(k1, k2)).read()
                    list = stdout.split()
                    if float(list[13]) < -9000:
                        print(stdout)

        # outFile.close()
