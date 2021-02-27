# rule split_overlap_chunks:
#     input:
#         fasta = "test/results/msa/"
#     output:
#         dir = directory("test/results/split_files/")
#     conda:
#         "env/pyfaidx.yaml"
#     script:
#         "../scripts/overlap_segs.py"



rule split_overlap_chunks:
    input:
        "test/results/msa/{sample}.msa.fasta"
    output:
        dir = directory("test/results/split_files/")
    conda:
        "env/pyfaidx.yaml"
    params:
        step = config["step"]
        overlap = config["overlap"]
    run:
        import os
        from pyfaidx import Fasta

        fasta_file = input[0]
        step = params.step
        overlap = params.overlap

        os.makedirs(output.dir)

        seqLength = len(fasta_file[0])

        remainder = seqLength % (step-overlap)
        chunkNumber = int(seqLength / (step-overlap))
        print("The sequences are splitted into {} chunks, and there are {} bp left".format())
