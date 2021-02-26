rule split_overlap_chunks:
    input:
        fasta = "test/mafft/sequences_aligned.fasta"
    output:
        dir = directory("test/splitFiles/")
    script:
        "scripts/overlap_segs.py"
