rule within_segments:
    input:
        datadir + "/{sample}/" + kmerCounter + "/forward_{seg}.kCountSorted"
    output:
        datadir + "/{sample}/primer3/forward_{seg}.calculated"
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

        output = open(output[0], "w")

        output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("kmer", "position", "kmerCount", "CG%", "Tm", "homodimer-dG", "hairpin-dG"))

        with open(input[0], "r") as input:
            for li in input:
                li = li.rstrip("\n")
                ls = li.split()
                kmer = ls[0]
                # primer3 tm calculation, but since RT alone does not have any
                # other substance than template and oligos. better use a simple
                # formule.
                # tm = primer3.calcTm(kmer,
                #                     mv_conc = params.monovalentConc,
                #                     dv_conc = params.divalentConc,
                #                     dntp_conc = params.dNTPConc,
                #                     dna_conc = params.dnaConc)
                kmer_CG = kmer.count("G") + kmer.count("C")
                kmer_AT = kmer.count("A") + kmer.count("T")
                tm = 64.9 + 41 * (kmer_CG-16.4)/(kmer_CG+kmer_AT)
                cgPercent = round(GC(kmer),3)
                homodimer_dg = primer3.calcHomodimer(kmer,
                                                     mv_conc = params.monovalentConc,
                                                     dv_conc = params.divalentConc,
                                                     dntp_conc = params.dNTPConc,
                                                     dna_conc = params.dnaConc).dg
                hairpin_dg = primer3.calcHairpin(kmer).dg
                output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(ls[0], wildcards.seg, ls[1], cgPercent, tm, homodimer_dg, hairpin_dg))
                # linguisticComplexity = 100 *

        output.close()

rule filterOut_within_segments:
    input:
        datadir + "/{sample}/primer3/forward_{seg}.calculated"
    output:
        datadir + "/{sample}/primer3/forward_{seg}.filtered"
    run:
        import pandas as pd
        import os

        df = pd.read_table(input[0], sep="\t", header=0)

        mean = float(df["Tm"].mean())
        std = float(df["Tm"].std())

        # MSSPE Deng et al. (2020), don't know the philosophy behind yet. with tm range from 60 - 70 no oligo could pass the criteria
        TmSeuilPlus = mean + 2*std
        TmSeuilLess = mean - 2*std

        df_filtered = df[(df["Tm"] >= TmSeuilLess) & (df["Tm"] <= TmSeuilPlus) & (df["CG%"] >= 40) & (df["CG%"] <= 60) & (df["hairpin-dG"] > -9000) & (df["homodimer-dG"] > -9000)]
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
