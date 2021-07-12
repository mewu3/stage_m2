rule deduplication:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        f"{dataDir}/{{sample}}/{{sample}}.uniq" if deduplication else []
    shell:
        "lib/cdhit/cd-hit-auxtools/cd-hit-dup -i {input} -o {output}"
