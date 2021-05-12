# viroExplorer
Pathogen viruses have an important impact in human history: epidemic and pandemic could cause great social and economic loss. Therefore it seems to us it is important to be able to follow their activities. 



This is a snakemake pipeline designed to construct oligonucleotide sequences for RT-PCR of group of study. Users are invited to modify the config.yaml file to choose the adequate parameters and define path for needed files such as: 
* nucl_gb.accession2taxid (ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/)
* ncbi_lineages_YYYY-MM-DD.csv (ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy, it is recommended to install via https://github.com/zyxue/ncbitax2lin) 
* etc. 



The philosophy of this pipeline is inspired by the work of Dr. Charles Chiu (https://github.com/chiulab/MSSPE-design), a primer design pipeline integrated in [1]. At this basic, the workflow was modified to adapt for environmental sample. 

basic: 
![image](https://user-images.githubusercontent.com/60400481/117315718-3d2bd480-ae88-11eb-9765-9503ff252b67.png)

clustering: same procedure than the basic, a step of clustering is added before multiple sequence alignment. 
![image](https://user-images.githubusercontent.com/60400481/117315764-474dd300-ae88-11eb-8637-b745fed86184.png)

variant: 
![image](https://user-images.githubusercontent.com/60400481/117317311-b11aac80-ae89-11eb-8527-976d7d91dcee.png)



Usage example (*It is recommended to create a conda environment for better tools management*): 

1. snakemake --cores 6 --snakefile snakefile-basic 
2. snakemake --cores 6 --snakefile snakefile-variant
3. snakemake --cores 6 --snakefile snakefile-clustering  (this one is not fully tested)

If you have few sequences (~ 3000 enterovirus genomes with a genome size around 7000 nt) you could use (1), but if the genome size is higher than 10 000 nt and the sequence number is greater than 5000 it is recommended to use (2). 
For more information about the snakemake API: https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html



[1] Deng, X., Achari, A., Federman, S. et al. Metagenomic sequencing with spiked primer enrichment for viral diagnostics and genomic surveillance. Nat Microbiol 5, 443–454 (2020). https://doi-org.inee.bib.cnrs.fr/10.1038/s41564-019-0637-9
