rule dskReverse:
    input:
        f"{dataDir}/{{sample}}/splitFiles/reverse{{seg}}.fasta"
    output:
        f"{dataDir}/{{sample}}/dsk/reverse{{seg}}.out"
    run: 
        import sys
        import os
        import re

        ### global variables 
        inputFile = input[0]
        outputFile = output[0]
        kmerCount = seqNumber = int(os.popen(f"egrep -c '^>' {inputFile}").read())
        filename = os.path.splitext(inputFile)[0]
        outputFile = open(outputFile, "w")

        ### functions 
        def calculate_kmer(input):

            dskH5 = f"{filename}.h5"
            os.system(
                f"lib/dsk/build/bin/dsk \
                -verbose 0 \
                -nb-cores 0 \
                -max-memory 5000 \
                -max-disk 0 \
                -out-compress 0 \
                -storage-type hdf5 \
                -kmer-size 13 \
                -abundance-min 2 \
                -abundance-max 2147483647 \
                -solidity-kind sum \
                -file {input} \
                -out {dskH5}"
            )

            dskh5txt = f"{filename}.kCount"
            os.system(
                f"lib/dsk/build/bin/dsk2ascii \
                -verbose 0 \
                -nb-cores 0 \
                -file {dskH5} -out {dskh5txt}"
            )

            dskh5txtSorted = f"{filename}.kCountSorted"
            os.system(
                f"sort -s -n -k 2 -nr {dskh5txt} > {dskh5txtSorted}"
            )

            dict = {}
            
            if os.stat(dskh5txtSorted).st_size != 0: 
                file = open(dskh5txtSorted, "r") 
                for line in file.readlines()[0:2]:
                    kmer = line.split()[0]
                    kmerCount = line.split()[1]
                    dict[kmer]=kmerCount
                file.close()

            return dict

        def parse_fasta(input): 

            dict={}

            file = open(input, "r")

            for line in file: 
                line = line.rstrip("\n")
                if line.startswith(">"): 
                    header = line 
                    dict[header] = ""
                else: 
                    dict[header] += line 

            file.close()

            return(dict)

        # ~ while kmerCount > seqNumber*0.1: 
           
            # ~ kmerDict = calculate_kmer(inputFile)
            # ~ kmerList = [key for key in kmerDict if kmerDict]     
            # ~ fastaDict = parse_fasta(inputFile)

            # ~ if kmerList: 
                # ~ for kmer in kmerList: 
                    # ~ boolean = [bool(re.search(kmer, fastaDict[key], re.I)) for key in fastaDict]
                    # ~ if True in boolean: 
                        # ~ kmerCount = int(kmerDict[kmer])
                        # ~ outputFile.write(f"{kmer}\t{kmerDict[kmer]}\n")
                        # ~ intermediate = open(f"{filename}.inter", "w")
                        # ~ for key in fastaDict: 
                            # ~ seq = fastaDict[key]
                            # ~ if re.search(kmer, seq, re.I) is None: 
                                # ~ intermediate.write(f"{key}\n{seq}\n")
                        # ~ intermediate.close()
                        # ~ break 
                    
            # ~ inputFile = f"{filename}.inter"
            
        while kmerCount > seqNumber*0.1: 
            
            kmerDict = calculate_kmer(inputFile)
            kmerList = [key for key in kmerDict if kmerDict] 
            fastaDict = parse_fasta(inputFile)

            if kmerList: 
                for kmer in kmerList: 
                    boolean = [bool(re.search(kmer, fastaDict[key], re.I)) for key in fastaDict]
                    if True in boolean: 
                        kmerCount = int(kmerDict[kmer])
                        outputFile.write(f"{kmer}\t{kmerDict[kmer]}\n")
                        intermediate = open(f"{filename}.inter", "w")
                        for key in fastaDict: 
                            seq = fastaDict[key]
                            if re.search(kmer, seq, re.I) is None: 
                                intermediate.write(f"{key}\n{seq}\n")
                        intermediate.close()
                        inputFile = f"{filename}.inter"
                        break 
                
        outputFile.close()
     
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.splitFiles.get(**wildcards).output[0]
    return expand(
        f"{dataDir}/{{sample}}/dsk/reverse{{seg}}.out",
        sample=wildcards.sample,
        seg=glob_wildcards(os.path.join(checkpoint_output, "reverse{seg}.fasta")).seg
    )

rule aggregate: 
	input: 
		aggregate_input
	output: 
		f"{dataDir}/{{sample}}/dsk/test.txt"
	shell: 
		"touch {output}"
