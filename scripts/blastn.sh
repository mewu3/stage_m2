#!/usr/bin/env bash

inputDir=$1

for fasta in $(ls $inputDir/reverse*.fasta)
do
    filename=$(basename $fasta)
    filename=${filename%.*}
    path=$(dirname ${fasta})
    blastn -task blastn-short \
    -db refseq_genomes -remote \
    -max_target_seqs 5 \
    -query $fasta \
    -out $path/${filename}.blastn.out \
    -outfmt '6 qacc sacc stitle pident length mismatch gapopen sstart send evalue'
done
