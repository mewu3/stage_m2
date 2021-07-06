rule clustering:
    input:
        dir = "/datater/wu/data/test_clustering/enterovirus"
    output:
        dir = "/datater/wu/data/test_clustering/drep_out"
    shell:
        "dRep cluster -p 6 -sa 0.95 -nc 0.85 -g {input} {output}"
