all = "/datater/wu/data/coronavirus/coronavirus.txt"

rule test:
    input:
        all
    output:
        "/datater/wu/data/genomeVariability/coronavirus"
    shell:
        "lib/FastANI/fastANI --rl {input[0]} --ql {input[0]} -o {output} --matrix -t 12"
