rule kmerFiltering1:
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

rule kmerFiltering2:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/allKmerCount.sorted.calculated.txt",
        input_kmerCounting1
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

        numSeq = int(os.popen(f"grep '^>' -c {input[1]}").read())
        numSeq_10 = int(numSeq * 0.1)

        df = pd.read_table(input[0], sep="\t", header=0, index_col=0)

        tmMean = float(df["Tm"].mean())
        tmStd = float(df["Tm"].std())
        TmSeuilPlus = tmMean + 2*tmStd
        TmSeuilLess = tmMean - 2*tmStd

        df_filtered = df[(df["Tm"] >= TmSeuilLess) & (df["Tm"] <= TmSeuilPlus) & (df["Tm"] <= TmMax) & (df["CG%"] >= GCUp) & (df["CG%"] <= GCDown) & (df["hairpin-dG"] > deltaG) & (df["homodimer-dG"] > deltaG) & (df["Entropy"] > LCSeuil) & (df["kmerCount"] > numSeq_10)]
        df_filtered = df_filtered.sort_values("kmerCount", ascending=False)
        df_filtered.to_csv(output[0], sep='\t', index=True)

rule kmerFiltering3:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/allKmerCount.sorted.calculated.filtered.txt"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.fasta"
    run:
        import os
        from Bio.Seq import Seq

        with open(output[0], "w") as fo:
            with open(input[0], "r") as fi:
                for l in fi.readlines()[1:]:
                    l = l.rstrip("\n").split()
                    id = l[0]
                    oligo = str(Seq(l[1]))
                    fo.write(f">{id}\n{oligo}\n")

rule kmerFiltering4:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.fasta"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.clustered.fasta"
    params:
        threads = config["thread"],
        memory = config["cd-hit"]["memory"]
    shell:
        "./lib/cdhit/cd-hit-est \
        -i {input} \
        -o {output} \
        -c 0.88 \
        -T {params.threads} \
        -M {params.memory}"

rule kmerFiltering5:
    input:
        file_refSeq
    output:
        multiext(
            f"{file_refSeq}.DB.", "1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"
        )
    params:
        out = f"{file_refSeq}.DB",
        threads = config["thread"]
    shell:
        """
        bowtie-build \
        --threads {params.threads} \
        {input[0]} {params.out}
        """

rule kmerFiltering6:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.clustered.fasta",
        lambda wildcards: expand(
            f"{file_refSeq}.DB.{{ext}}",
            ext = ["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"]
        )
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.clustered.bowtie"
    params:
        refDB = f"{file_refSeq}.DB",
        threads = config["thread"]
    shell:
        """
        bowtie \
        -x {params.refDB} -f {input[0]} \
        -v 0 -l 7 -a \
        --sam \
        -p {params.threads} \
        {output}
        """

rule kmerFiltering7:
    input:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.clustered.bowtie",
        f"{dataDir}/{{sample}}/{{kmerSize}}/allKmerCount.sorted.calculated.filtered.txt"
    output:
        f"{dataDir}/{{sample}}/{{kmerSize}}/intermediate/allKmerCount.sorted.calculated.filtered.clustered.spec.fasta"
    params:
        kmerSize = lambda wildcards: config["kmerSize"][wildcards.kmerSize]
    run:
        import os
        import pandas as pd

        specific = os.popen(f"samtools view -f 4 {input[0]} | cut -f1").read().split("\n") # -f unmapped -F matched
        specific = list(filter(None, specific))
        specific = list(map(lambda x: int(x), specific))

        df = pd.read_table(input[1], sep="\t", header=0, index_col=0)
        df = df[df.index.isin(specific)]

        with open(output[0], "w") as fo:
            for i, r in df.iterrows():
                seq = r["oligo"]
                fo.write(f">{i}\n{seq}\n")
