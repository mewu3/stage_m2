library("tidyverse")
library("gridExtra")

### load input files
filePath = "MSSPE-basic"
virus = "enterovirusCurated" 
kmerSize = "15" 
datadir <- sprintf("/datater/wu/data/%s/%s/kmer%s/evaluation", filePath, virus, kmerSize)

kmerCount <- list.files(path=datadir, pattern = "tsv", full.names = TRUE)
# specieCoverage = paste(datadir, "allOligo.set.coverage", sep="/")
specieCoverage = paste(datadir, "allKmerCount.sorted.calculated.filtered.allOligo.set.coverage", sep="/")

### global variables
# totolSeq = 3265 # enterovirus
# totolSeq = 75959 # coronavirus
totolSeq = 119 # enterovirusCurated
# totolSeq = 258 # coronavirusCurated

### Basic & Clustering ######################################################

### plot oligonucleotides count coverage
# for (file in kmerCount){
#     filename <- tools::file_path_sans_ext(file)
#     filename <- paste(filename, ".png", sep="")
#     table <- read_delim(file=file, delim="\t", col_names=TRUE)
#     table[,"kmerCoverage"] <- table[,"kmerCount"]/totolSeq*100
#     table$position <- paste(table$start, table$end, sep="-")
#     table$position <- factor(table$position, levels=table$position)
#     plot <- ggplot(table, aes(x=position, y=kmerCoverage)) +
#       geom_col() +
#       theme(axis.text.x = element_text(angle=90, size=10)) +
#       ylim(0,100) +
#       xlab("position")
#       ylab("kmer Coverage")
#     ggsave(plot, file=filename, width=24, height=12, units="cm", scale=2)
# }

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
ggsave(file=filename, arrangeGrob(plot1, plot2), width=24, height=12, units="cm", scale=2)

# per species
# species = unique(table$species)
# for (specie in species){
#     specieTable = table %>% filter(species == specie)
#     filename <- tools::file_path_sans_ext(specieCoverage)
#     filename <- paste(filename, specie, sep=".")
#     filename <- paste(filename, ".png", sep="")
#     plot1 <- ggplot(specieTable, aes(x=position, y=speciesCoverage, fill=species)) +
#         geom_col() +
#         theme(axis.text.x = element_text(angle=90, size=10)) +
#         ylim(0,100) +
#         labs(x='position', y='Frequency of species', title='acutal species coverage')
#     plot2 <- ggplot(specieTable, aes(x=position, y=allSpeciesCoverage, fill=species)) +
#         geom_col() +
#         theme(axis.text.x = element_text(angle=90, size=10)) +
#         ylim(0,100) +
#         labs(x='position', y='Frequency of species', title='all species coverage')
#     ggsave(arrangeGrob(plot1, plot2), file=filename, width=24, height=12, units="cm", scale=2)
# }

### variant ##################################################################

# ### plot oligonucleotides count coverage
# for (file in kmerCount){
#     filename <- tools::file_path_sans_ext(file)
#     filename <- paste(filename, ".png", sep="")
#     table <- read_delim(file=file, delim="\t", col_names=TRUE)
#     table[,"kmerCount"] <- table[,"kmerCount"]/totolSeq*100
#     table$position <- paste(table$`chunk-start`, table$`chunk-end`, sep="-")
#     table$position <- factor(table$position, levels=table$position)
#     plot <- ggplot(table, aes(x=position, y=kmerCount)) +
#       geom_col() +
#       theme(axis.text.x = element_text(angle=90, size=10)) +
#       ylim(0,100) +
#       xlab("position")
#       ylab("kmer Coverage")
#     ggsave(plot, file=filename, width=24, height=12, units="cm", scale=2)
# }

# ### plot oligonucleotides species coverage
# table = read_delim(file=specieCoverage, delim="\t", col_names=FALSE)
# names(table) = c("start", "end", "chunk-start", "chunk-end", "species", "speciesCount", "totalCount")
# table[, "speciesCoverage"] = table[,"speciesCount"]/totolSeq*100
# table[, "allSpeciesCoverage"] = table[,"totalCount"]/totolSeq*100
# table <- arrange(table, start)
# table[, "position"] <- paste(table$`chunk-start`, table$`chunk-end`, sep="-")
# table$position <- factor(table$position, levels=unique(table$position))

# # all together
# plot1 <- ggplot(table, aes(x=position, y=speciesCoverage, fill=species)) +
#     geom_col() +
#     theme(axis.text.x = element_text(angle=90, size=10)) +
#     labs(x='position', y='Frequency of species', title='acutal species coverage')
# plot2 <- ggplot(table, aes(x=position, y=allSpeciesCoverage, fill=species)) +
#     geom_col() +
#     theme(axis.text.x = element_text(angle=90, size=10)) +
#     labs(x='position', y='Frequency of species', title='all species coverage')
# filename <- tools::file_path_sans_ext(specieCoverage)
# filename <- paste(filename, "allSpecies.png", sep=".")
# ggsave(file=filename, arrangeGrob(plot1, plot2), width=24, height=12, units="cm", scale=2)

# # per species
# species = unique(table$species)
# for (specie in species){
#     specieTable = table %>% filter(species == specie)
#     filename <- tools::file_path_sans_ext(specieCoverage)
#     filename <- paste(filename, specie, sep=".")
#     filename <- paste(filename, ".png", sep="")
#     plot1 <- ggplot(specieTable, aes(x=position, y=speciesCoverage, fill=species)) +
#         geom_col() +
#         theme(axis.text.x = element_text(angle=90, size=10)) +
#         labs(x='position', y='Frequency of species', title='acutal species coverage')
#     plot2 <- ggplot(specieTable, aes(x=position, y=allSpeciesCoverage, fill=species)) +
#         geom_col() +
#         theme(axis.text.x = element_text(angle=90, size=10)) +
#         labs(x='position', y='Frequency of species', title='all species coverage')
#     ggsave(arrangeGrob(plot1, plot2), file=filename, width=24, height=12, units="cm", scale=2)
# }