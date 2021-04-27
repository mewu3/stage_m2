rule specific1:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allKmerCount.sorted.calculated.filtered.txt"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.fasta"
    run:
        import os

        outFile = open(output[0], "w")
        input1Open = open(input[0], "r")
        for li in input1Open.readlines()[1:]:
            li = li.rstrip("\n")
            ls = li.split()
            id = ls[0]
            oligo = ls[1]
            outFile.write(f">p{id}\n{oligo}\n")
        input1Open.close()
        outFile.close()

# rule blastDB:
#     input:
#         refSeq
#     output:
#         multiext(
#             f"{refSeqFile}.DB.", "nhr", "nin", "nog", "nsd", "nsi", "nsq"
#         )
#     params:
#         out = f"{refSeqFile}.DB"
#     shell:
#         "makeblastdb -in {input} -out {params.out} \
#         -dbtype nucl \
#         -parse_seqids"
#
# rule blastn_short_prok:
#     input:
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.fasta"
#     output:
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.blastn"
#     params:
#         refDB = f"{refSeqFile}.DB",
#         maxTargetSeq = config["blast"]["max-target-seqs"],
#         thread = config["blast"]["thread"]
#     shell:
#         """
#         blastn -task blastn-short \
#         -db {params.refDB} \
#         -word_size 7 \
#         -evalue 1000 \
#         -dust no \
#         -max_target_seqs {params.maxTargetSeq} \
#         -num_threads {params.thread} \
#         -query {input} \
#         -out {output} \
#         -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore'
#         """
#         # 19/04 evalue dust and word size parameters fallows guide: https://www-ncbi-nlm-nih-gov.inee.bib.cnrs.fr/blast/BLAST_guide.pdf
#         # -outfmt '6 qacc sacc stitle pident length mismatch gapopen sstart send evalue'
#
# rule filterOutUnSpecifickmer:
#     input:
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.blastn",
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/allKmerCount.sorted.calculated.filtered.txt"
#     output:
#         f"{dataDir}/{{sample}}/kmer{kmerSize}/allKmerCount.sorted.calculated.filtered.spec.txt"
#     params:
#         kmerLength = config["jellyfish"]["kmer-size"]
#     run:
#         from Bio import SeqIO
#         import pandas as pd
#
#         kmerSize = int(params.kmerLength)
#
#         df_blast = pd.read_table(input[0], sep="\t")
#
#         df_blast_filtered = df_blast[(df_blast.iloc[:,3] == kmerSize) & (df_blast.iloc[:,4] == 0)] # 20/04
#         # 19/04 df_blast_filtered = df_blast[(df_blast.iloc[:,3] > kmerSize *0.8) & (df_blast.iloc[:,4] > kmerSize * 0.2)]
#
#         nonSpecificId = df_blast_filtered.iloc[:,0].tolist()
#         nonSpecificId = set(nonSpecificId)
#         nonSpecificId = [int(x.lstrip("p")) for x in nonSpecificId]
#
#         df = pd.read_table(input[1], sep="\t", header=0, index_col=0)
#         df_filtered = df[~df.index.isin(nonSpecificId)]
#         df_filtered.to_csv(output[0], sep='\t', index=True)

rule specific2:
    input:
        refSeq
    output:
        multiext(
            f"{refSeqFile}.DB.", "1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"
        )
    params:
        out = f"{refSeqFile}.DB",
        threads = 12
    shell:
        """
        bowtie-build \
        --threads {params.threads} \
        {input} {params.out}
        """

rule specific3:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.fasta",
        lambda wildcards: expand(
            f"{refSeqFile}.DB.{{ext}}",
            ext = ["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"]
        )
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.bowtie"
    params:
        refDB = f"{refSeqFile}.DB",
        threads = 12
    shell:
        """
        bowtie \
        -x {params.refDB} -f {input[0]} \
        -v 0 -l 7 -a \
        --sam \
        -p {params.threads} \
        {output}
        """

rule specific4:
    input:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.bowtie",
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allKmerCount.sorted.calculated.filtered.txt"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/allKmerCount.sorted.calculated.filtered.spec.txt"
    params:
        kmerLength = config["jellyfish"]["kmer-size"]
    run:
        import os
        import pandas as pd

        specific = os.popen(f"samtools view -f 4 {input[0]} | cut -f1").read().split("\n") # -f unmapped -F matched
        specific = list(map(lambda x: x.lstrip("p"), specific))
        specific = list(filter(None, specific))
        specific = list(map(lambda x: int(x), specific))

        df = pd.read_table(input[1], sep="\t", header=0, index_col=0)
        df_filtered = df[df.index.isin(specific)]
        df_filtered.to_csv(output[0], sep='\t', index=True)
