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


# ~ checkpoint splitFilteredTable:
    # ~ input:
        # ~ f"{dataDir}/{{sample}}/filtering/allOligos_reverse.filtered"
    # ~ output:
        # ~ directory(f"{dataDir}/{{sample}}/checkSpecifity")
        # ~ dynamic(f"{dataDir}/{{sample}}/checkSpecifity/reverse{{seg}}_filtered.tsv")
    # ~ run:
        # ~ import os
        # ~ from Bio.Seq import Seq
        # ~ import pandas as pd

        # ~ df = pd.read_table(input[0], sep="\t", header=0)
        # ~ for seg, table in df.groupby("position"):
            # ~ outputFile = open(f"{dataDir}/{wildcards.sample}/checkSpecifity/reverse{seg}_filtered.tsv", "w")
            # ~ # outputFile = open(f"{output[0]}/reverse{seg}_filtered.tsv", "w")
            # ~ table.to_csv(outputFile, sep="\t", index=False)
            # ~ outputFile.close()

# ~ rule splitTableToFasta:
    # ~ input:
        # ~ f"{dataDir}/{{sample}}/checkSpecifity/reverse{{seg}}_filtered.tsv"
    # ~ output:
        # ~ f"{dataDir}/{{sample}}/checkSpecifity/reverse{{seg}}_filtered.fasta"
    # ~ run:
        # ~ import os

        # ~ outFile = open(output[0], "w")
        # ~ n = 0
        # ~ with open(input[0], "r") as file:
            # ~ for li in file.readlines()[1:21]:
                # ~ li = li.rstrip("\n")
                # ~ ls = li.split()
                # ~ oligo = ls[7]
                # ~ outFile.write(f">p{n}|position {ls[0]}|count {ls[1]}|CG% {ls[2]}|Tm {ls[3]}|homodimer_dG {ls[4]}|hairpin_dG {ls[5]}|LC {ls[6]}\n{oligo}\n")
                # ~ n += 1
        # ~ outFile.close()

# ~ def aggregate_input(wildcards):
    # ~ checkpoint_out = checkpoints.splitFilteredTable.get(**wildcards).output[0]
    # ~ return expand(
        # ~ f"{dataDir}/{{sample}}/checkSpecifity/reverse{{seg}}_filtered.fasta",
        # ~ sample = wildcards.sample,
        # ~ seg = glob_wildcards(os.path.join(checkpoint_out, "reverse{seg}_filtered.tsv")).seg
    # ~ )

# ~ rule aggregateAllToFasta:
	# ~ input:
		# ~ aggregate_input
	# ~ output:
		# ~ f"{dataDir}/{{sample}}/checkSpecifity/reverse{{seg}}_filtered.fasta"

# rule RNAtoDNA:
#     input:
#         "{dataDir}/{{sample}}/{{sample}}.uniq".format(dataDir=dataDir)filtering
#     output:
#         "{dataDir}/{{sample}}/{{sample}}.uniq.dna".format(dataDir=dataDir)
#     shell:
#         "python3 scripts/RNAtoDNA.py {input} {output}"


### vmath alignment: examination needed for the output in the future
# rule mktree_index:
#     input:
#         "{dataDir}/{{sample}}/{{sample}}.uniq.dna".format(dataDir=dataDir)
#     output:filtering
#         multiext(
#             "{dataDir}/{{sample}}/{{sample}}.index.".format(dataDir=dataDir), "al1", "bck", "bwt", "des", "lcp", "llv", "ois", "prj", "sds", "skp", "ssp", "sti1", "suf", "tis"
#         )
#     params:
#         index = "{dataDir}/{{sample}}/{{sample}}.index".format(dataDir=dataDir)
#     shell :
#         "lib/vmatch-2.3.1/mkvtree -db {input} -dna -pl -allout -indexname {params.index}"
#
# rule vmatch:
#     input:
#         "{dataDir}/{{sample}}/checkSpecifity/allOligos_reverse.fasta".format(dataDir=dataDir)
#     output:
#         "{dataDir}/{{sample}}/checkSpecifity/vmatch.out".format(dataDir=dataDir)
#     params:
#         index = "{dataDir}/{{sample}}/{{sample}}.index".format(dataDir=dataDir)
#     shell:filtering
#         "lib/vmatch-2.3.1/vmatch -complete -p -v -showdesc 10 -s 120 abbrev -noevalue -noscore -noidentity -q {input} {params.index} > {output}"

### last alignment: output produce nothing
# rule lastdb_index:
#     input:
#         "{dataDir}/{{sample}}/{{sample}}.uniq.dna".format(dataDir=dataDir)
#     output:
#         multiext(
#             "{dataDir}/{{sample}}/{{sample}}.index.".format(dataDir=dataDir), "bck", "des", "prj", "sds", "ssp", "suf", "tis"
#         )
#     params:
#         index = "{dataDir}/{{sample}}/{{sample}}.index".format(dataDir=dataDir)
#     shell:
#         "lib/last-1179/src/lastdb {params.index} {input}"
#
# rule lastal_alignment:
#     input:
#         "{dataDir}/{{sample}}/checkSpecifity/allOligos_reverse.fasta".format(dataDir=dataDir)
#     output:
#         "{dataDir}/{{sample}}/checkSpecifity/lastal.out".format(dataDir=dataDir)
#     params:
#         index = "{dataDir}/{{sample}}/{{sample}}.index".format(dataDir=dataDir)
#     shell:
#         "lib/last-1179/src/lastal {params.index} {input} > {output}"

### bbmap: don't have position information
# rule bbmap:
#     input:
#         query = "{dataDir}/{{sample}}/checkSpecifity/allOligos_reverse.fasta".format(dataDir=dataDir),
#         ref = "{dataDir}/{{sample}}/{{sample}}.uniq.dna".format(dataDir=dataDir)
#     output:
#         "{dataDir}/{{sample}}/checkSpecifity/seal.out".format(dataDir=dataDir),
#         "{dataDir}/{{sample}}/checkSpecifity/seal.out.stats".format(dataDir=dataDir)
#     params:
#         refPath = "{dataDir}/{{sample}}/".format(dataDir=dataDir)
#     shell:
#         "lib/bbmap/seal.sh in={input.query} ref={input.ref} out={output[0]} stats={output[1]} k=6 rcomp=t mm=f"

### good old blast blatn short alignement, if remote is used don't need to build blastDB
# rule blastDB:
#     input:
#         refSeq
#     output:
#         multiext(
#             f"{dataDir}/{{sample}}/checkSpecifity{kmerSize}/references.DB.", "nhr", "nin", "nog", "nsd", "nsi", "nsq"
#         )
#     params:
#         out = f"{dataDir}/{{sample}}/checkSpecifity{kmerSize}/references.DB"
#     shell:
#         "makeblastdb -in {input} -out {params.out} \
#         -dbtype nucl \
#         -parse_seqids"
#
# rule blastn_short_prok:
#     input:
#         f"{dataDir}/{{sample}}/checkSpecifity{kmerSize}/allOligos_reverse.filtered.fasta"
#     output:
#         f"{dataDir}/{{sample}}/checkSpecifity{kmerSize}/allOligos_reverse.filtered.blastn.out"
#     params:
#         refDB = f"{dataDir}/{{sample}}/checkSpecifity{kmerSize}/references.DB",
#         evalue = config["blast"]["evalue"],
#         wordSize = config["blast"]["word-size"],
#         maxTargetSeq = config["blast"]["max-target-seqs"]
#     shell:
#         "blastn -task blastn-short \
#         -db {params.refDB} \
#         -word_size {params.wordSize} \
#         -max_target_seqs {params.maxTargetSeq} \
#         -query {input} \
#         -out {output} \
#         -outfmt '6 qacc sacc stitle pident length mismatch gapopen sstart send evalue'"
#
# rule filterOutUnSpecifickmer:
#     input:
#         f"{dataDir}/{{sample}}/checkSpecifity{kmerSize}/allOligos_reverse.filtered.blastn.out",
#         f"{dataDir}/{{sample}}/checkSpecifity{kmerSize}/allOligos_reverse.filtered.fasta"
#     output:
#         f"{dataDir}/{{sample}}/checkSpecifity{kmerSize}/allOligos_reverse.filtered.spec.fasta"
#     run:
#         from Bio import SeqIO
#
#         fastaFile = input[1]
#         blastnFile = input[0]
#         outFasta = output[0]
#
#         blastnFileOpen = open(blastnFile, "r")
#         fastaFileOpen = open(fastaFile, "r")
#         outputOpen = open(outFasta, "w")
#
#         pIDSet = set()
#         for line in blastnFileOpen:
#             line = line.rstrip("\n")
#             list = line.split("\t")
#             pID = list[0]
#             identity = float(list[3])
#             if identity > 50:
#                 pIDSet.add(pID)
#
#         seqID = []
#         for line in fastaFileOpen:
#             line = line.rstrip("\n")
#             if line.startswith(">"):
#                 header = line.split("|")[0].lstrip(">").strip()
#                 if header not in pIDSet:
#                     seqID.append(header)
#
#         record_dict = SeqIO.to_dict(SeqIO.parse(fastaFile, "fasta"))
#         for ID in seqID:
#             fasta = record_dict[ID].format("fasta")
#             outputOpen.write(f"{fasta}")
#
#         blastnFileOpen.close()
#         fastaFileOpen.close()
#         outputOpen.close()

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
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.fasta"
    output:
        f"{dataDir}/{{sample}}/kmer{kmerSize}/intermediate/allKmerCount.sorted.calculated.filtered.bowtie"
    params:
        refDB = f"{refSeqFile}.DB",
        threads = 12
    shell:
        """
        bowtie \
        -x {params.refDB} -f {input} \
        -v 0 -l 7 \
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
