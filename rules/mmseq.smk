rule mmseq2: # cluster input sequences
    input:
        "test/enterovirus_50.fasta"
    output:
        "test/enterovirus_50_clsutered.fasta"
    log:
        "test/mmseq2.log"
    shell:
        "mmseqs createdb {input} | \
         mmseq2 easy-linclust {output} "
