library("tidyverse")
library("gridExtra")

### load input files
# datadir <- "/home/meijun/Documents/server/data/MSSPE-basic/enterovirus/kmer13/evaluation"
# datadir <- "/datater/wu/data/MSSPE-basic/enterovirus/kmer13/evaluation"
# datadir <- "/home/meijun/Documents/server/data/MSSPE-basic/enterovirus/kmer15/evaluation"
datadir <- "/datater/wu/data/MSSPE-clustering/enterovirus/kmer15/evaluation"
kmerCount <- list.files(path=datadir, pattern = "tsv", full.names = TRUE)
specieCoverage = paste(datadir, "allOligo.set.coverage", sep="/")

### global variables
totolSeq = 3265

### plot oligonucleotides count coverage
for (file in kmerCount){
    filename <- tools::file_path_sans_ext(file)
    filename <- paste(filename, ".png", sep="")
    table <- read_delim(file=file, delim="\t", col_names=TRUE)
    table[,"kmerCoverage"] <- table[,"kmerCount"]/totolSeq*100
    table$position <- paste(table$start, table$end, sep="-")
    table$position <- factor(table$position, levels=table$position)
    plot <- ggplot(table, aes(x=position, y=kmerCoverage)) +
      geom_col() +
      theme(axis.text.x = element_text(angle=90, size=10)) +
      ylim(0,100) +
      xlab("position")
      ylab("kmer Coverage")
    ggsave(plot, file=filename)
}

### plot oligonucleotides species coverage
table = read_delim(file=specieCoverage, delim="\t", col_names=FALSE)
names(table) = c("start", "end", "species", "speciesCount", "totalCount")
table[, "speciesCoverage"] = table[,"speciesCount"]/totolSeq*100
table[, "allSpeciesCoverage"] = table[,"totalCount"]/totolSeq*100
table <- arrange(table, start)
table$position <- paste(table$start, table$end, sep="-")
table$position <- factor(table$position, levels=unique(table$position))

# all together
plot1 <- ggplot(table, aes(x=position, y=speciesCoverage, fill=species)) +
    geom_col() +
    theme(axis.text.x = element_text(angle=90, size=10)) +
    ylim(0,100) +
    labs(x='position', y='Frequency of species', title='acutal species coverage')
plot2 <- ggplot(table, aes(x=position, y=allSpeciesCoverage, fill=species)) +
    geom_col() +
    theme(axis.text.x = element_text(angle=90, size=10)) +
    ylim(0,100) +
    labs(x='position', y='Frequency of species', title='all species coverage')
filename <- tools::file_path_sans_ext(specieCoverage)
filename <- paste(filename, "allSpecies.png", sep=".")
ggsave(file=filename, arrangeGrob(plot1, plot2))

# per species
species = unique(table$species)
for (specie in species){
    specieTable = table %>% filter(species == specie)
    filename <- tools::file_path_sans_ext(specieCoverage)
    filename <- paste(filename, specie, sep=".")
    filename <- paste(filename, ".png", sep="")
    plot1 <- ggplot(specieTable, aes(x=position, y=speciesCoverage, fill=species)) +
        geom_col() +
        theme(axis.text.x = element_text(angle=90, size=10)) +
        ylim(0,100) +
        labs(x='position', y='Frequency of species', title='acutal species coverage')
    plot2 <- ggplot(specieTable, aes(x=position, y=allSpeciesCoverage, fill=species)) +
        geom_col() +
        theme(axis.text.x = element_text(angle=90, size=10)) +
        ylim(0,100) +
        labs(x='position', y='Frequency of species', title='all species coverage')
    ggsave(arrangeGrob(plot1, plot2), file=filename)
}
