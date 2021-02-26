rule kmc:
    input:
        dir = "test/splitFiles"
    output:
        dir = directory("test/kmerCounting/kmc3")
    # log:
    #     "test/kmerCounting/kmc3/log.txt"
    conda:
        "envs/kmc3.yaml"
    params:
        kmerSize = config["kmc3"]["kmersize"]
    script:
        "scripts/kmc3.py"
