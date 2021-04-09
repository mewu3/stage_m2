rule calculateKmer_DSK_h5:
    input:
        f"{dataDir}/{{sample}}/splitFiles/reverse{{seg}}.fasta"
    output:
        f"{dataDir}/{{sample}}/kmerCounting/reverse{{seg}}.jf"
    shell:
        "jellyfish count -m 13 -s 100M -o {output} {input} "

rule calculateKmer_DSK_txt:
    input:
        f"{dataDir}/{{sample}}/kmerCounting/reverse{{seg}}.jf"
    output:
        f"{dataDir}/{{sample}}/kmerCounting/reverse{{seg}}.kCount"
    params:
        nbCores = config["dsk2ascii"]["nb-cores"],
    shell:
        "jellyfish dump -c {input} > {output}"

rule calculateKmer_DSK_txt_sort:
    input:
        f"{dataDir}/{{sample}}/kmerCounting/reverse{{seg}}.kCount"
    output:
        f"{dataDir}/{{sample}}/kmerCounting/reverse{{seg}}.kCountSorted"
    run:
        with open(input[0], "r") as input:
            rows = input.readlines()
            sorted_rows = sorted(rows, key = lambda x: int(x.split()[1]), reverse=True)
            with open(output[0], "w") as output:
                for row in sorted_rows:
                        output.write(row)
