# viroExplorer
Pathogen viruses have an important impact in human history: epidemic and pandemic could cause great social and economic loss. Therefore it seems to us it is important to be able to follow their activities. 



This is a snakemake pipeline designed to construct oligonucleotide sequences for RT-PCR of group of study. Users are invited to modify the config.yaml file to choose the adequate parameters and define path for needed files such as: 
* nucl_gb.accession2taxid (ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/)
* ncbi_lineages_YYYY-MM-DD.csv (ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy, it is recommended to install via https://github.com/zyxue/ncbitax2lin) 
* etc. 



The philosophy of this pipeline is inspired by the work of Dr. Charles Chiu (https://github.com/chiulab/MSSPE-design), a primer design pipeline integrated in [1]. At this basic, the workflow was modified to adapt for environmental sample. 

basic: 
![image](https://user-images.githubusercontent.com/60400481/117315718-3d2bd480-ae88-11eb-9765-9503ff252b67.png)

basic + clustering: same procedure than the basic, a step of clustering is added before multiple sequence alignment. clustering can be enabled via config.yaml file. 
![image](https://user-images.githubusercontent.com/60400481/117315764-474dd300-ae88-11eb-8637-b745fed86184.png)

variant: 
![image](https://user-images.githubusercontent.com/60400481/117317311-b11aac80-ae89-11eb-8527-976d7d91dcee.png)



Usage example (*It is recommended to create a conda environment for better tools management*): 

1. snakemake --cores 6 --snakefile snakefile-basic 
2. snakemake --cores 6 --snakefile snakefile-variant

With ~ 3000 enterovirus complet genomes (genome size around 7000 nt), basic approch could be done within 30 mins. If the NL (N sequence number, L sequence length) is important, it is recommended to use (2) ou enable the clustering option. 
For more information about the snakemake API: https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html



[1] Deng, X., Achari, A., Federman, S. et al. Metagenomic sequencing with spiked primer enrichment for viral diagnostics and genomic surveillance. Nat Microbiol 5, 443–454 (2020). https://doi-org.inee.bib.cnrs.fr/10.1038/s41564-019-0637-9
