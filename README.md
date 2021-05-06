# viroExplorer
Pathogen virus have a important impact in human history: epidemic and pandemic could cause great social and economic loss. Therefor it seems to us it's important to be able to follow their activities. 

This is a snakemake pipeline designed to construct oligonucleotide sequences for RT-PCR of study group. Users are invited to modify the config.yaml file to choose the adequate parameters and define path for needed files such as: 
* nucl_gb.accession2taxid (ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/)
* ncbi_lineages_YYYY-MM-DD.csv (ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy, it's recommended to use https://github.com/zyxue/ncbitax2lin) 
* etc. 

The pilosophie of this pipeliine is inspiered by the work of Dr. Charles Chiu (https://github.com/chiulab/MSSPE-design), a primer design pipeliine integrated in [1]. At this basic, the workflow was modified to adapt for environmental sample. 



[1] Deng, X., Achari, A., Federman, S. et al. Metagenomic sequencing with spiked primer enrichment for viral diagnostics and genomic surveillance. Nat Microbiol 5, 443–454 (2020). https://doi-org.inee.bib.cnrs.fr/10.1038/s41564-019-0637-9
