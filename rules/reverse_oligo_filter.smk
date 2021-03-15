rule reverse_parameters:
    input:
        datadir + "/{sample}/" + kmerCounter + "/reverse{seg}.kCountSorted"
    output:
        datadir + "/{sample}/filtering/reverse{seg}.calculated"
    params:
        monovalentConc = config["oligotm"]["monovalent-conc"],
        divalentConc = config["oligotm"]["divalent-conc"],
        dNTPConc = config["oligotm"]["dNTP-conc"],
        dnaConc = config["oligotm"]["dna-conc"],
        thermodynamicPara = config["oligotm"]["thermodynamic-para"],
        saltCorrelation = config["oligotm"]["salt-correlation"],
    run:
        import os
        from Bio.SeqUtils import GC
        import primer3
        from Bio.Seq import Seq

        output = open(output[0], "w")

        output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("oligo", "position", "kmerCount", "CG%", "Tm", "homodimer-dG", "hairpin-dG"))

        with open(input[0], "r") as input:
            for li in input:
                li = li.rstrip("\n")
                ls = li.split()
                kmer = Seq(ls[0])
                oligo = str(kmer.reverse_complement())
                # primer3 tm calculation, but since RT alone does not have any
                # other substance than template and oligos. better use a simple
                # formule.
                # tm = primer3.calcTm(kmer,
                #                     mv_conc = params.monovalentConc,
                #                     dv_conc = params.divalentConc,
                #                     dntp_conc = params.dNTPConc,
                #                     dna_conc = params.dnaConc)
                kmer_CG = oligo.count("G") + oligo.count("C")
                kmer_AT = oligo.count("A") + oligo.count("T")
                tm = 64.9 + 41 * (kmer_CG-16.4)/(kmer_CG+kmer_AT)
                cgPercent = round(GC(oligo),3)
                homodimer_dg = primer3.calcHomodimer(oligo,
                                                     mv_conc = params.monovalentConc,
                                                     dv_conc = params.divalentConc,
                                                     dntp_conc = params.dNTPConc,
                                                     dna_conc = params.dnaConc).dg
                hairpin_dg = primer3.calcHairpin(oligo).dg
                output.write("{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\n".format(oligo, wildcards.seg.split("-")[1], ls[1], cgPercent, tm, homodimer_dg, hairpin_dg))

        output.close()


### calculate complexity with nessie, input format require fasta
rule reverse_toFasta:
    input:
        datadir + "/{sample}/filtering/reverse{seg}.calculated"
    output:
        datadir + "/{sample}/filtering/reverse{seg}.fa"
    run:
        import os

        outFile = open(output[0], "w")

        n = 0
        with open(input[0], "r") as file:
            next(file)
            for li in file:
                li = li.rstrip("\n")
                ls = li.split()
                outFile.write(">p{}|position {}|count {}|CG% {}|Tm {}|homodimer_dG {}|hairpin_dG {}\n{}\n".format(n, ls[1], ls[2], ls[3], ls[4], ls[5], ls[6], ls[0]))
                n += 1

        outFile.close()

rule reverse_linguistic_complexity:
    input:
        datadir + "/{sample}/filtering/reverse{seg}.fa"
    output:
        datadir + "/{sample}/filtering/reverse{seg}.nessieOut"
    shell:
        "lib/nessie/nessie -I {input} -O {output} -L"

rule reverse_add_LC:
    input:
        datadir + "/{sample}/filtering/reverse{seg}.nessieOut",
        datadir + "/{sample}/filtering/reverse{seg}.fa"
    output:
        datadir + "/{sample}/filtering/reverse{seg}.calculated2"
    run:
        input2 = open(input[1], "r")
        outputFile = open(output[0], "w")

        dict_id_lc = {}

        with open(input[0], "r") as input1 :
            next(input1)
            for li1 in input1:
                li1 = li1.rstrip("\n")
                if li1.startswith(">"):
                    id = li1.split("|")[0].lstrip(">")
                    dict_id_lc[id] = 0
                else:
                    LC = li1.split(":")[1]
                    dict_id_lc[id] = LC

        ### write to fasta format
        # for li2 in input2:
        #     li2 = li2.rstrip("\n")
        #     seq = ""
        #     newHeader = ""
        #     if li2.startswith(">"):
        #         header2 = li2
        #         id_ = li2.split("|")[0].lstrip(">")
        #         for id in dict_id_lc:
        #             if id_ == id :
        #                 newHeader = "{} |LC {:.3f}".format(header2, dict_id_lc[id])
        #                 outputFile.write(newHeader+"\n")
        #     else:
        #         seq = li2
        #         outputFile.write(seq + "\n")

        ### write to tsv format
        outputFile.write("position\tkmerCount\tCG%\tTm\thomodimer-dG\thairpin-dG\tLC\toligo\n")
        for li2 in input2:
            li2 = li2.rstrip("\n")
            if li2.startswith(">"):
                ls = li2.split("|")
                id_ = ls[0].lstrip(">")
                position = ls[1].split()[1]
                count = ls[2].split()[1]
                cg = ls[3].split()[1]
                tm = ls[4].split()[1]
                homodimer_dG = ls[5].split()[1]
                hairpin_dG = ls[6].split()[1]
                outputFile.write("{}\t{}\t{}\t{}\t{}\t{}\t".format(position, count, cg, tm, homodimer_dG, hairpin_dG))
                for id in dict_id_lc:
                    if id_ == id :
                        lc = float(dict_id_lc[id_])
                        outputFile.write("{:.3f}\t".format(lc))
            else:
                seq = li2
                outputFile.write("{}\n".format(seq))

        outputFile.close()

### error message, merge in interactiv commande line
# def aggregate_reverseInput(wildcards):
#     outputdir = datadir + "/{sample}/filtering/"
#     return expand(datadir + "/{sample}/filtering/reverse{seg}.calculated2",
#                   sample = wildcards.sample,
#                   seg = glob_wildcards(os.path.join(outputdir, "reverse{seg}.calculated2")).seg)
#
# rule aggregate_allReverseOligo:
#     input:
#         aggregate_reverseInput
#     output:
#         datadir + "/{sample}/filtering/allOligos_reverse.calculated2"
#     shell:
#         """
#         echo -e "position\tkmerCount\tCG%\tTm\thomodimer-dG\thairpin-dG\tLC\tkmer" > {output}
#         """
#         """
#         sed -s 1d {input} >> {output}
#         """

rule reverse_filtering:
    input:
        datadir + "/{sample}/filtering/allOligos_reverse.calculated2"
    output:
        datadir + "/{sample}/filtering/allOligos_reverse.filtered"
    run:
        import pandas as pd
        import os

        df = pd.read_table(input[0], sep="\t", header=0)

        mean = float(df["Tm"].mean())
        std = float(df["Tm"].std())

        # MSSPE Deng et al. (2020), don't know the philosophy behind yet. with tm range from 60 - 70 no oligo could pass the criteria
        TmSeuilPlus = mean + 2*std
        TmSeuilLess = mean - 2*std

        df_filtered = df[(df["Tm"] >= TmSeuilLess) & (df["Tm"] <= TmSeuilPlus) & (df["CG%"] >= 40) & (df["CG%"] <= 60) & (df["hairpin-dG"] > -9000) & (df["homodimer-dG"] > -9000) & (df["LC"] >= 0.75)]
        df_filtered = df_filtered.sort_values("position", ascending=False)
        df_filtered.to_csv(output[0], sep='\t', index=False)


# def aggregate_input(wildcards):
#     checkpoint_output = datadir + "/{sample}/primer3_Tm/"
#     return expand(datadir + "/{sample}/primer3_Tm/forward_{seg}.TmFilterered",
#                   sample = wildcards.sample,
#                   seg = glob_wildcards(os.path.join(checkpoint_output, "forward_{seg}.TmFilterered")).seg)

# rule merge_all_kmers:
#     input:
#         aggregate_input
#     output:
#         datadir + "/{sample}/primer3_Tm/all_oligo"
#     shell:
#         "sed -s 1d {input} >> {output}"


#
# rule dimer_mfeprimer:
#     input:
#         datadir + "/{sample}/primer3_Tm/all_oligo.fasta"
#     output:
#         datadir + "/{sample}/MFEprimer/all_oligo_hairpin_dimer"
#     params:
#         deltaG = config["mfeprimer"]["deltaG-cutoff"],
#         Tm = config["mfeprimer"]["Tm-cutoff"],
#         score = config["mfeprimer"]["score-cutoff"],
#         monovalent = config["mfeprimer"]["monovalent-conc"],
#         divalent = config["mfeprimer"]["divalent-conc"],
#         NPT = config["mfeprimer"]["dNTP-conc"],
#         oligo = config["mfeprimer"]["oligo-conc"]
#     shell:
#         "lib/mfeprimer-3.2.1-linux-amd64 dimer -i {input} -o {output} \
#         -t {params.Tm} \
#         -d {params.deltaG} \
#         -s {params.score} \
#         --mono {params.monovalent} \
#         --diva {params.divalent} \
#         --dntp {params.NPT} \
#         --oligo {params.oligo}"

#
# rule dimer_ntthal:
#     input:
#         datadir + "/{sample}/primer3_ntthal/allOligo_hairpin"
#     output:
#         datadir + "/{sample}/primer3_ntthal/allOligo_hairpin_dimer"
#     run:
#         import numpy as np
#         import os
#         # import pandas as pd
#
#         inputFile = open(input[0], "r")
#         outputFile = open(output[0], "w")
#
#         lines = []
#
#         for li in inputFile:
#             li = li.rstrip("\n")
#             ls = li.split()
#             lines.append(ls)
#
#         for k1 in np.array(lines).T[0]:
#             for k2 in np.array(lines).T[0]:
#                 stdout = os.popen("lib/primer3/src/ntthal -s1 {} -s2 {}".format(k1, k2)).read()
#                 ls = stdout.split()
#                 deltaG = ls[13]
#                 if float(deltaG) < -9000:
#                     # print(k1, k2, deltaG)
#                     outputFile.write("{}\t{}\t{}\n".format(k1, k2, deltaG))
#                     # outputFile.write(stdout)
#                 # outputFile.write(stdout)
#
#         outputFile.close()
