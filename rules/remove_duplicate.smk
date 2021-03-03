rule seqkit: # remove duplicate sequences ######################################
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        # "{}/duplicate_removed/{sample}.uniq.fasta".format(datadir)
        datadir + "/duplicate_removed/{sample}.uniq.fasta"
    conda:
        "envs/seqkit.yaml"
    log:
        # "{}/duplicate_removed/seqkit.{sample}.log".format(datadir)
        datadir + "/duplicate_removed/{sample}.uniq.fasta"
    shell:
        "cat {input} | seqkit rmdup -s -o {output}"
