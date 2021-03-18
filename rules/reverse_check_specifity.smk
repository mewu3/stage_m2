rule allToFasta:
    input:
        "{dataDir}/{{sample}}/filtering/allOligos_reverse.filtered".format(dataDir=dataDir)
    output:
        "{dataDir}/{{sample}}/checkSpecifity/allOligos_reverse.fasta".format(dataDir=dataDir)
    run:
        import os
        from Bio.Seq import Seq

        outFile = open(output[0], "w")

        n = 0
        with open(input[0], "r") as file:
            next(file)
            for li in file:
                li = li.rstrip("\n")
                ls = li.split()
                # oligo = ls[7]
                oligo = str(Seq(ls[7]).reverse_complement())
                outFile.write(">p{}|position {}|count {}|CG% {}|Tm {}|homodimer_dG {}|hairpin_dG {}|LC {}\n{}\n".format(n, ls[0], ls[1], ls[2], ls[3], ls[4], ls[5], ls[6], oligo))
                n += 1

        outFile.close()

rule RNAtoDNA:
    input:
        "{dataDir}/{{sample}}/{{sample}}.uniq".format(dataDir=dataDir)
    output:
        "{dataDir}/{{sample}}/{{sample}}.uniq.dna".format(dataDir=dataDir)
    shell:
        "python3 scripts/RNAtoDNA.py {input} {output}"


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

### good old blast blatn short liangment
rule blastDB:
    input:
        "{dataDir}/{{sample}}/{{sample}}.uniq".format(dataDir=dataDir)
    output:
        multiext(
            "{dataDir}/{{sample}}/{{sample}}.uniq.db.".format(dataDir=dataDir), "nhr", "nin", "nog", "nsd", "nsi", "nsq"
        )
    params:
        title = "Making_blastDB",
        out = "{dataDir}/{{sample}}/{{sample}}.uniq.db".format(dataDir=dataDir)
    shell:
        "makeblastdb -in {input} -out {params.out} \
        -parse_seqids \
        -title {params.title} \
        -dbtype nucl"

rule blastn_short:
    input:
        "{dataDir}/{{sample}}/checkSpecifity/allOligos_reverse.fasta".format(dataDir=dataDir)
    output:
        "{dataDir}/{{sample}}/checkSpecifity/blastn.out".format(dataDir=dataDir)
    params:
        ref = "{dataDir}/{{sample}}/{{sample}}.uniq.db".format(dataDir=dataDir)
    shell:
        "blastn -task blastn-short \
        -query {input} \
        -db {params.ref} \
        -out {output} \
        -outfmt 7 "

# rule parse_blastOutput:
#     input:
#         "{dataDir}/{{sample}}/checkSpecifity/blastn.out".format(dataDir=dataDir)
#     output:
#         "{dataDir}/{{sample}}/checkSpecifity/blastn.out.parsed".format(dataDir=dataDir)
#     run:
#         from collections import defaultdict
#
#         dict_count = defaultdict(lambda: defaultdict(int))
#         dict_posit = defaultdict(lambda: defaultdict(list))
#
#         with open(input[0], "r") as inputFile:
#             for line in inputFile:
#                 line = line.rstrip("\n")
#                 if not line.startswith("#"):
#                     ls = line.split()
#                     qId = ls[0]
#                     refId = ls[1]
#                     position = [ls[-4], ls[-3]]
#                     dict_count[qId][refId] += 1
#                     dict_posit[qId][refId].append(position)
#
#         for qId in dict_count:
#             print(dict_count[qId])
#             # for refId in dict_count[qId]:
#                 # print(dict_count[qId][refId])


rule blastn_short_nr:
    input:
        "{dataDir}/{{sample}}/checkSpecifity/allOligos_reverse.fasta".format(dataDir=dataDir)
    output:
        "{dataDir}/{{sample}}/checkSpecifity/blastn_nt.out".format(dataDir=dataDir)
    shell:
        "blastn -task blastn-short \
        -query {input} \
        -db nt \
        -out {output} \
        -outfmt 7 \
        -remote"
