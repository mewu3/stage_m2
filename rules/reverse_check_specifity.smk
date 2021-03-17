rule allToFasta:
    input:
        "{dataDir}/{{sample}}/filtering/allOligos_reverse.filtered".format(dataDir=dataDir)
    output:
        "{dataDir}/{{sample}}/filtering/allOligos_reverse.fasta".format(dataDir=dataDir)
    run:
        import os
        from Bio.Seq import Seq

        outFile = open(output[0], "w")

        n = 0
        with open(input[0], "r") as file:
            next(file)
            for li in file:
                li = li.rstrip("\n")
                ls = li.split()
                # oligo = ls[7]
                oligo = str(Seq(ls[7]).reverse_complement())
                outFile.write(">p{}|position {}|count {}|CG% {}|Tm {}|homodimer_dG {}|hairpin_dG {}|LC {}\n{}\n".format(n, ls[0], ls[1], ls[2], ls[3], ls[4], ls[5], ls[6], oligo))
                n += 1

        outFile.close()
