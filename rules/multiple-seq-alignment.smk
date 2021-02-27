if config["mafft"]["fmodel"] == "on":
    config["mafft"]["fmodel"] = "--fmodel"
else:
    config["mafft"]["fmodel"] = ""

if config["mafft"]["clustalout"] == "on":
    config["mafft"]["clustalout"] = "--clustalout"
else:
    config["mafft"]["clustalout"] = ""

if config["mafft"]["inputorder"] == "on":
    config["mafft"]["inputorder"] = "--inputorder"
else :
    config["mafft"]["inputorder"] = ""

if config["mafft"]["reorder"] == "on":
    config["mafft"]["reorder"] = "--reorder"
else:
    config["mafft"]["reorder"] = ""

if config["mafft"]["treeout"] == "on":
    config["mafft"]["treeout"] = "--treeout"
else :
    config["mafft"]["treeout"] = ""

if config["mafft"]["quiet"] == "on":
    config["mafft"]["quiet"] = "--quiet"
else:
    config["mafft"]["quiet"] = ""

rule mafft:
    input:
        "test/results/duplicate_removed/{sample}.uniq.fasta"
    output:
        "test/results/msa/{sample}.msa.fasta"
    log:
        "test/results/msa/mafft.{sample}.log"
    conda:
        "envs/mafft.yaml"
    params:
        threads = config["mafft"]["threads"],
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
        "mafft \
        --thread {params.threads} \
        {params.algorithm} \
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
