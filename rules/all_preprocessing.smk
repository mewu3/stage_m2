rule seqkit: # remove duplicate sequences ######################################
    input:
        sampleDir
    output:
        "{}/duplicate_removed/{}.uniq".format(dataDir, sample)
    conda:
        "envs/seqkit.yaml"
    shell:
        "cat {input} | seqkit rmdup -s -o {output}"

rule mafft: # multiple sequences alignment #####################################
    input:
        "{}/duplicate_removed/{}.uniq".format(dataDir, sample)
    output:
        "{}/msa/{}.msa".format(dataDir, sample)
    log:
        "{}/msa/mafft.{}.log".format(dataDir, sample)
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

checkpoint splitFiles:
    input:
        "{}/msa/{}.msa".format(dataDir, sample)
    output:
        directory("{}/{}/splitFiles".format(dataDir, sample))
    params:
        step = config["step"],
        overlap = config["overlap"]
    run:
        #!/usr/bin/env python3.8
        from pyfaidx import Fasta
        import os

        inputFile = input[0]
        os.makedirs(output[0])
        step = params.step
        overlap = params.overlap

        fasta = Fasta(inputFile)
        seqLength = len(fasta[0])

        remainder = seqLength % (step-overlap)
        chunkNumber = int(seqLength / (step-overlap))
        print("The sequences are splitted into "+str(chunkNumber)+" chunks, and there are "+str(remainder)+" bp left.")

        if remainder <= step/2: # primux fasta_tile_overlap.pl
            newStep = int(remainder/chunkNumber) + 1 + step
            print("Changing step size from {} to {} so there will be no remainder.".format(step, newStep))

        chunks = [[i,i+newStep] for i in range(0, seqLength, newStep-overlap)]
        chunks[-1][1] = len(fasta[0])

        for chunk in chunks:

            seg = str(chunk[0])+"-"+str(chunk[1])

            f1 = open(output[0] + "/forward{}.fasta".format(seg), "w")
            f2 = open(output[0] + "/reverse{}.fasta".format(seg), "w")

            for id in fasta.keys() :
                segment = fasta[id][chunk[0]:chunk[1]]
                forward = str(segment[:50])
                reverse = str(segment[-50:])
                f1.write("> {} |{}-{} \n {} \n".format(fasta[id].long_name, str(chunk[0]), str(chunk[1]), forward))
                f2.write("> {} |{}-{} \n {} \n".format(fasta[id].long_name, str(chunk[0]), str(chunk[1]), reverse))

            f1.close()
            f2.close()