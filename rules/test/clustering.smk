rule clustering:
    input:
        f"{dataDir}/{{sample}}/{{sample}}.uniq" if deduplication else lambda wildcards: config["samples"][wildcards.sample]
    output:
        f"{dataDir}/{{sample}}/{{sample}}.uniq.cluster" if deduplication else f"{dataDir}/{{sample}}/{{sample}}.cluster"
    params:
        threads = config["thread"],
        memory = config["cd-hit"]["memory"]
    shell:
        "./lib/cdhit/cd-hit-est \
        -i {input} \
        -o {output} \
        -c 0.95 -aS 0.85 -aL 0.85 \
        -T {params.threads} \
        -M {params.memory}"
