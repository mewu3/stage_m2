# rule cdhit: # cluster input sequences
#     input:
#         "{samples_path}/{samples}.fasta"
#     output:
#         dir = directory("{results_dir}/cdhit")
#     log:
#         "{results_dir}/cdhit.log"
#     shell:
#         "lib/cdhit/cdhit -i {input} -o {samples} "

rule seqkit:
    input:
        "test/{sample}.fasta"
    output:
        "test/results/duplicate_removed/{sample}.uniq.fasta"
    conda:
        "envs/seqkit.yaml"
    log:
        "test/results/duplicate_removed/seqkit.{sample}.log"
    shell:
        "cat {input} |seqkit rmdup -s -o {output}"
