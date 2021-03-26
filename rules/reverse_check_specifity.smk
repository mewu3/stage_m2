rule splitFilteredTable:
    input:
        f"{dataDir}/{{sample}}/filtering/allOligos_reverse.filtered"
    output:
        # directory(f"{dataDir}/{{sample}}/checkSpecifity")
        dynamic(f"{dataDir}/{{sample}}/checkSpecifity/reverse{{seg}}_filtered.tsv")
    run:
        import os
        from Bio.Seq import Seq
        import pandas as pd

        df = pd.read_table(input[0], sep="\t", header=0)
        for seg, table in df.groupby("position"):
            outputFile = open(f"{dataDir}/{wildcards.sample}/checkSpecifity/reverse{seg}_filtered.tsv", "w")
            # outputFile = open(f"{output[0]}/reverse{seg}_filtered.tsv", "w")
            table.to_csv(outputFile, sep="\t", index=False)
            outputFile.close()

rule splitTableToFasta:
    input:
        f"{dataDir}/{{sample}}/checkSpecifity/reverse{{seg}}_filtered.tsv"
    output:
        f"{dataDir}/{{sample}}/checkSpecifity/reverse{{seg}}_filtered.fasta"
    run:
        import os

        outFile = open(output[0], "w")
        n = 0
        with open(input[0], "r") as file:
            next(file)
            for li in file:
                li = li.rstrip("\n")
                ls = li.split()
                oligo = ls[7]
                outFile.write(f">p{n}|position {ls[0]}|count {ls[1]}|CG% {ls[2]}|Tm {ls[3]}|homodimer_dG {ls[4]}|hairpin_dG {ls[5]}|LC {ls[6]}\n{oligo}\n")
                n += 1
        outFile.close()

# def aggregate_input(wildcards):
#     checkpoint_out = checkpoints.splitFilteredTable.get(**wildcards).output[0]
#     return expand(
#         f"{dataDir}/{{sample}}/checkSpecifity/reverse{{seg}}_filtered.fasta",
#         sample = wildcards.sample,
#         seg = glob_wildcards(os.path.join(checkpoint_out, "reverse{seg}_filtered.tsv")).seg
#     )

# rule RNAtoDNA:
#     input:
#         "{dataDir}/{{sample}}/{{sample}}.uniq".format(dataDir=dataDir)
#     output:
#         "{dataDir}/{{sample}}/{{sample}}.uniq.dna".format(dataDir=dataDir)
#     shell:
#         "python3 scripts/RNAtoDNA.py {input} {output}"


### vmath alignment: examination needed for the output in the future
# rule mktree_index:
#     input:
#         "{dataDir}/{{sample}}/{{sample}}.uniq.dna".format(dataDir=dataDir)
#     output:
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
#     shell:
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
#         "{dataDir}/{{sample}}/checkSpecifity/references.fasta".format(dataDir=dataDir)
#     output:
#         multiext(
#             "{dataDir}/{{sample}}/checkSpecifity/references.DB.".format(dataDir=dataDir), "nhr", "nin", "nog", "nsd", "nsi", "nsq"
#         )
#     params:
#         title = "Building_blastDB",
#         out = "{dataDir}/{{sample}}/checkSpecifity/references.DB".format(dataDir=dataDir)
#     shell:
#         "makeblastdb -in {input} -out {params.out} \
#         -parse_seqids \
#         -title {params.title} \
#         -dbtype nucl"

rule blastn_short:
    input:
        f"{dataDir}/{{sample}}/checkSpecifity/reverse{{seg}}_filtered.fasta"
    output:
        f"{dataDir}/{{sample}}/checkSpecifity/reverse{{seg}}_filtered.blatn.out"
    shell:
        "blastn -task blastn-short \
        -db refseq_genomes -remote \
        -query {input} \
        -out {output} \
        -outfmt '6 qacc sacc stitle pident length mismatch gapopen sstart send evalue'"
