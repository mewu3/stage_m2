rule seqkit: # remove duplicate sequences ######################################
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        datadir + "/duplicate_removed/{sample}.uniq"
    conda:
        "envs/seqkit.yaml"
    shell:
        "cat {input} | seqkit rmdup -s -o {output}"
