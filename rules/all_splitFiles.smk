rule splitFiles: # split MSA fasta file into seperates files #########
    input:
        datadir + "/msa/{sample}.msa"
    output:
        # dynamic(datadir + "/{sample}/splitFiles/forward{seg}.fasta"),
        dynamic(datadir + "/{sample}/splitFiles/reverse{seg}.fasta")
    params:
        step = config["step"],
        overlap = config["overlap"]
    run:
        #!/usr/bin/env python3.8
        from pyfaidx import Fasta

        inputFile = input[0]
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

            f1 = open(datadir + "/" + wildcards.sample + "/splitFiles/forward{}.fasta".format(seg), "w")
            f2 = open(datadir + "/" + wildcards.sample + "/splitFiles/reverse{}.fasta".format(seg), "w")

            for id in fasta.keys() :
                segment = fasta[id][chunk[0]:chunk[1]]
                forward = str(segment[:50])
                reverse = str(segment[-50:])
                f1.write("> {} |{}-{} \n {} \n".format(fasta[id].long_name, str(chunk[0]), str(chunk[1]), forward))
                f2.write("> {} |{}-{} \n {} \n".format(fasta[id].long_name, str(chunk[0]), str(chunk[1]), reverse))

            f1.close()
            f2.close()
