rule dskReverse:
    input:
        f"{dataDir}/{{sample}}/splitFiles/reverse{{seg}}.fasta"
    output:
        f"{dataDir}/{{sample}}/dsk/reverse{{seg}}.out"
    shell: 
        "python3 scripts/kmerSelect_iterac.py {input} {output}"

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
