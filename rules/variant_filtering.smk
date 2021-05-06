rule filtering1:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.txt"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/allKmerCount.sorted.calculated.txt"
    params:
        monovalentConc = config["filter-oligotm"]["monovalent-conc"],
        divalentConc = config["filter-oligotm"]["divalent-conc"],
        dNTPConc = config["filter-oligotm"]["dNTP-conc"],
        dnaConc = config["filter-oligotm"]["dna-conc"],
        thermodynamicPara = config["filter-oligotm"]["thermodynamic-para"],
        saltCorrelation = config["filter-oligotm"]["salt-correlation"],
    run:
        import os
        from Bio.SeqUtils import GC
        import primer3
        from Bio.Seq import Seq
        import collections
        import math

        def estimate_shannon_entropy(dna_sequence):

            m = len(dna_sequence)
            bases = collections.Counter([tmp_base for tmp_base in dna_sequence])

            shannon_entropy_value = 0
            for base in bases:
                n_i = bases[base]
                p_i = n_i / float(m)
                entropy_i = p_i * (math.log(p_i, 2))
                shannon_entropy_value += entropy_i

            return shannon_entropy_value * -0.5 # nessie shannon entropy

        # from: https://github.com/bioinform/somaticseq/blob/e6f5b1c6b98b324d418
        # 407154392778164215a65/somaticseq/utilities/linguistic_sequence_complex
        # ity.py
        # According to: https://doi.org/10.1093/bioinformatics/18.5.679
        # Assume 4 different nucleotides

        def max_vocabularies(seq_length):

            counts = 0
            k = 1
            while k <= seq_length:

                if 4**k < (seq_length - k + 1):
                    counts = counts + 4**k
                else:
                    counts = counts + (seq_length-k+1 + 1) * (seq_length-k+1 - 1 + 1)/2
                    break

                k += 1

            return counts

        def LC(sequence):

            sequence = sequence.upper()

            if not 'N' in sequence:

                number_of_subseqs     = 0
                seq_length            = len(sequence)
                max_number_of_subseqs = max_vocabularies(seq_length)

                for i in range(1, seq_length+1):

                    #max_vocab_1 = 4**i
                    #max_vocab_2 = seq_length - i + 1
                    set_of_seq_n = set()

                    for n, nth_base in enumerate(sequence):

                        if n+i <= len(sequence):
                            sub_seq = sequence[n:n+i]
                            set_of_seq_n.add( sub_seq )

                    num_uniq_subseqs  = len(set_of_seq_n)
                    number_of_subseqs = number_of_subseqs + num_uniq_subseqs

                lc = number_of_subseqs/max_number_of_subseqs

            else:
                lc = float('nan')

            return lc

        output1Open = open(output[0], "w")
        output1Open.write(f"\toligo\tkmerCount\tCG%\tEntropy\tTm\thomodimer-dG\thairpin-dG\n")

        input1Open = open(input[0], "r")
        n=0
        for li in input1Open:
            li = li.rstrip("\n")
            ls = li.split()
            oligo = Seq(ls[0])
            oligo = str(oligo.reverse_complement())
            # primer3 tm calculation, but since RT alone does not have any
            # other substance than template and oligos. better use a simple
            # formule.
            count = ls[1]
            cgPercent = round(GC(oligo),3)
            entropy = LC(oligo)
            kmer_CG = oligo.count("G") + oligo.count("C")
            kmer_AT = oligo.count("A") + oligo.count("T")
            tm = 64.9 + 41 * (kmer_CG-16.4)/(kmer_CG+kmer_AT)
            homodimer_dg = primer3.calcHomodimer(oligo,
                                                 mv_conc = params.monovalentConc,
                                                 dv_conc = params.divalentConc,
                                                 dntp_conc = params.dNTPConc,
                                                 dna_conc = params.dnaConc).dg
            hairpin_dg = primer3.calcHairpin(oligo).dg
            entropy = estimate_shannon_entropy(oligo)
            output1Open.write(f"{n}\t{oligo}\t{count}\t{cgPercent:.3f}\t{entropy:.3f}\t{tm:.3f}\t{homodimer_dg:.3f}\t{hairpin_dg:.3f}\n")
            n += 1

        input1Open.close()
        output1Open.close()

rule filtering2:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/allKmerCount.sorted.calculated.txt"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/allKmerCount.sorted.calculated.filtered.txt"
    params:
        deltaG = config["homodimer-deltaG"],
        GCUp = config["GC-upper"],
        GCDown = config["GC-lower"],
        TmMax = config["Tm-max"],
        LCSeuil = config["linguistic-complexity"]
    run:
        import pandas as pd
        import os

        deltaG = params.deltaG
        GCUp = params.GCUp
        GCDown = params.GCDown
        TmMax = params.TmMax
        LCSeuil = params.LCSeuil

        df = pd.read_table(input[0], sep="\t", header=0, index_col=0)

        # entropyThreshold = float(df["Entropy"].quantile(0.75))
        # entropyThreshold = float(df["Entropy"].median())
        # entropyThreshold = float(df["Entropy"].mean())
        tmMean = float(df["Tm"].mean())
        tmStd = float(df["Tm"].std())

        # MSSPE Deng et al. (2020), don't know the logic behind yet. with tm
        # range from 60 - 70 no oligo could pass the criteria
        TmSeuilPlus = tmMean + 2*tmStd
        TmSeuilLess = tmMean - 2*tmStd

        df_filtered = df[(df["Tm"] >= TmSeuilLess) & (df["Tm"] <= TmSeuilPlus) & (df["Tm"] <= TmMax) & (df["CG%"] >= GCUp) & (df["CG%"] <= GCDown) & (df["hairpin-dG"] > deltaG) & (df["homodimer-dG"] > deltaG) & (df["Entropy"] > LCSeuil)]
        df_filtered = df_filtered.sort_values("kmerCount", ascending=False)
        df_filtered.to_csv(output[0], sep='\t', index=True)
