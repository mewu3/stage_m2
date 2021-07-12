def input(wildcards): 
    if deduplication and clustering: 
        return expand(
            f"{dataDir}/{{sample}}/{{sample}}.uniq.cluster{clusterIdentity}",
            sample = wildcards.sample
        )
    elif deduplication: 
        return expand(
            "{dataDir}/{{sample}}/{{sample}}.uniq",
            sample = wildcards.sample
        ) 
    elif clustering: 
        return expand(
            f"{dataDir}/{{sample}}/{{sample}}.cluster{clusterIdentity}",
            sample = wildcards.sample
        )
    else: 
        return config["samples"][wildcards.sample]


rule MSA:
    input:
        input
    output:
        f"{dataDir}/{{sample}}/{{sample}}.msa"
    log:
        f"{dataDir}/{{sample}}/log/mafft.log"
    conda:
        "envs/mafft.yaml"
    params:
        threads = config["thread"],
        algorithm = config["mafft"]["algorithm"],
        op = config["mafft"]["op"],
        ep = config["mafft"]["ep"],
        bl = config["mafft"]["bl"],
        jtt = config["mafft"]["jtt"],
        tm = config["mafft"]["tm"],
        fmodel = config["mafft"]["fmodel"],
        clustalout = config["mafft"]["clustalout"],
        inputorder = config["mafft"]["inputorder"],
        reorder = config["mafft"]["reorder"],
        treeout = config["mafft"]["treeout"],
        quiet = config["mafft"]["quiet"]
    shell:
        "mafft {params.algorithm} \
        --thread {params.threads} \
        --op {params.op} \
        --ep {params.ep} \
        --bl {params.bl} \
        --jtt {params.jtt} \
        --tm {params.tm} \
        {params.fmodel} \
        {params.clustalout} \
        {params.inputorder} \
        {params.reorder} \
        {params.treeout} \
        {params.quiet} \
        {input} > {output} \
        2> {log}"
