rule removeDuplicateSeq:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        f"{dataDir}/{{sample}}/{{sample}}.uniq"
    shell:
        "lib/cdhit/cd-hit-auxtools/cd-hit-dup -i {input} -o {output}"

if clustering :
    rule clustering:
        input:
            f"{dataDir}/{{sample}}/{{sample}}.uniq"
        output:
            f"{dataDir}/{{sample}}/{{sample}}.uniq.cluster{clusterIdentity}"
        params:
            identity = config["cd-hit"]["identity"],
            threads = config["thread"],
            memory = config["cd-hit"]["memory"]
        shell:
            "./lib/cdhit/cd-hit-est \
            -i {input} \
            -o {output} \
            -c {params.identity} \
            -T {params.threads} \
            -M {params.memory}"

    rule MSA:
        input:
            f"{dataDir}/{{sample}}/{{sample}}.uniq.cluster{clusterIdentity}"
        output:
            f"{dataDir}/{{sample}}/{{sample}}.uniq.cluster{clusterIdentity}.msa"
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
else:
    rule MSA:
        input:
            lambda wildcards: config["samples"][wildcards.sample] if config["curated"] else f"{dataDir}/{{sample}}/{{sample}}.uniq"
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
            """
            mafft {params.algorithm} \
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
            2> {log}
            """

checkpoint splitIntoOverlappingWindows:
    input:
        f"{dataDir}/{{sample}}/{{sample}}.uniq.cluster{clusterIdentity}.msa" if clustering else f"{dataDir}/{{sample}}/{{sample}}.msa"
    output:
        directory(f"{dataDir}/{{sample}}/splitFiles")
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
            step = newStep

        chunks = [[i,i+step] for i in range(0, seqLength, step-overlap)]
        chunkNumber=len(chunks)
        print(f"Final chunk number: {chunkNumber}")

        for chunk in chunks:

            if chunk[-1] > seqLength:
                chunk[-1] = seqLength

            seg = str(chunk[0])+"-"+str(chunk[1])

            # f1 = open(output[0] + f"/forward{seg}.fasta", "w")
            f2 = open(output[0] + f"/reverse{seg}.fasta", "w")

            for id in fasta.keys() :
                segment = fasta[id][chunk[0]:chunk[1]]
                # forward = str(segment[:50])
                reverse = str(segment[-50:])
                # f1.write("> {}|{}-{}\n{}\n".format(fasta[id].long_name, str(chunk[0]), str(chunk[1]), forward))
                f2.write(f">{fasta[id].long_name}\n{reverse}\n")

            # f1.close()
            f2.close()
