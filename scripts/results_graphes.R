library("tidyverse")
library("gridExtra")
library("reshape2")

### load input files
datadir = "/home/meijun/Documents/server/data/MSSPE-basic"
virus = c("enterovirus", "coronavirus", "coronavirusWithoutSarsCov2", "enterovirusCurated", "coronavirusCurated", "coronavirusCuratedWithoutSarsCov2")
kmerSize = c("13", "15")
combinations = expand.grid(virus, kmerSize)

for (row in 1:nrow(combinations)){
  virus = toString(combinations[row,1])
  kmerSize = toString(combinations[row,2])
  workdir <- sprintf("%s/%s/kmer%s/evaluation", datadir, virus, kmerSize)
  if (dir.exists(workdir)){
    ### kmerCoverage counting ###
    oligoCount_before <- sprintf("%s/allOligo_before.tsv", workdir)
    oligoCount_before <- read.table(file=oligoCount_before, 
                                    sep="\t", 
                                    header=TRUE)
    oligoCount_before <-oligoCount_before %>% select(start, end, kmerCount)
    oligoCount_after1 <- sprintf("%s/allOligo_after1.tsv", workdir)
    oligoCount_after1 <- read.table(file=oligoCount_after1, 
                                    sep="\t", 
                                    header=TRUE)
    oligoCount_after1 <- oligoCount_after1 %>% select(start, end, kmerCount)
    oligoCount_after2 <- sprintf("%s/allOligo_after2.tsv", workdir)
    oligoCount_after2 <- read.table(file=oligoCount_after2, 
                                    sep="\t", 
                                    header=TRUE)
    oligoCount_after2 <- oligoCount_after2 %>% select(start, end, kmerCount)
    oligoCount_after3 <- sprintf("%s/allOligo_after3.tsv", workdir)
    oligoCount_after3 <- read.table(file=oligoCount_after3, 
                                    sep="\t", 
                                    header=TRUE)
    oligoCount_after3 <-oligoCount_after3 %>% select(start, end, kmerCount)
    oligoCount <- left_join(x=oligoCount_before, y=oligoCount_after1, by=c("start", "end"), suffix=c(".a0", ".a1"))
    oligoCount <- left_join(x=oligoCount, y=oligoCount_after2, by=c("start", "end"), suffix=c(".a1", ".a2"))
    oligoCount <- left_join(x=oligoCount, y=oligoCount_after3, by=c("start", "end"), suffix=c(".a2", ".a3"))
    oligoCount[is.na(oligoCount)] <- 0
    oligoCount <- oligoCount %>% arrange(start, decreasing=TRUE)
    oligoCount$position <- paste(oligoCount$start, oligoCount$end, sep="-")
    oligoCount$position <- factor(oligoCount$position, levels=oligoCount$position)
    oligoCount <- oligoCount[, !(names(oligoCount) %in% c("start", "end"))]
    dfm <- pivot_longer(oligoCount, -position, names_to="variable", values_to="value")
    plot <- ggplot(dfm, aes(x=position,y=value)) +
      geom_bar(aes(fill=variable), stat = "identity", position = "dodge") + 
      theme(axis.text.x = element_text(angle=90, size=5)) +
      ylim(0,100) +
      xlab("position") + 
      ylab("kmer Coverage")
    print(plot)
    
    ### species coverage ###
    finalSpeciesCoverage = sprintf("%s/allOligo.set.coverage", workdir)
    finalSpeciesCoverage = read.table(file=finalSpeciesCoverage, 
                                      sep="\t", 
                                      header=FALSE)
    names(finalSpeciesCoverage) = c("start", "end", "species", "speciesCount", "totalCount")
    
  }
}

datadir <- sprintf("%s/%s/%s/kmer%s/evaluation", datadir, filePath, virus, kmerSize)

kmerCount <- list.files(path=datadir, pattern = "tsv", full.names = TRUE)
specieCoverage = paste(datadir, "allOligo.set.coverage", sep="/")
# specieCoverage = paste(datadir, "allKmerCount.sorted.calculated.filtered.allOligo.set.coverage", sep="/")

### global variables
# totolSeq = 3265 # enterovirus
# totolSeq = 75959 # coronavirus
# totolSeq = 119 # enterovirusCurated
totolSeq = 258 # coronavirusCurated

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
# table = read_delim(file=specieCoverage, delim="\t", col_names=FALSE)
# names(table) = c("start", "end", "species", "speciesCount", "totalCount")
# table[, "speciesCoverage"] = table[,"speciesCount"]/totolSeq*100
# table[, "allSpeciesCoverage"] = table[,"totalCount"]/totolSeq*100
# table <- arrange(table, start)
# table$position <- paste(table$start, table$end, sep="-")
# table$position <- factor(table$position, levels=unique(table$position))
#
# # all together
# plot1 <- ggplot(table, aes(x=position, y=speciesCoverage, fill=species)) +
#     geom_col() +
#     theme(axis.text.x = element_text(angle=90, size=10)) +
#     ylim(0,100) +
#     labs(x='position', y='Frequency of species', title='acutal species coverage')
# plot2 <- ggplot(table, aes(x=position, y=allSpeciesCoverage, fill=species)) +
#     geom_col() +
#     theme(axis.text.x = element_text(angle=90, size=10)) +
#     ylim(0,100) +
#     labs(x='position', y='Frequency of species', title='all species coverage')
# filename <- tools::file_path_sans_ext(specieCoverage)
# filename <- paste(filename, "allSpecies.png", sep=".")
# ggsave(file=filename, arrangeGrob(plot1, plot2), width=24, height=12, units="cm", scale=2)

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

### plot oligonucleotides species coverage
table = read_delim(file=specieCoverage, delim="\t", col_names=FALSE)
names(table) = c("start", "end", "chunk-start", "chunk-end", "species", "speciesCount", "totalCount")
table[, "speciesCoverage"] = table[,"speciesCount"]/totolSeq*100
table[, "allSpeciesCoverage"] = table[,"totalCount"]/totolSeq*100
table <- arrange(table, start)
table[, "position"] <- paste(table$`chunk-start`, table$`chunk-end`, sep="-")
table$position <- factor(table$position, levels=unique(table$position))

# all together
plot1 <- ggplot(table, aes(x=position, y=speciesCoverage, fill=species)) +
    geom_col() +
    theme(axis.text.x = element_text(angle=90, size=10)) +
    labs(x='position', y='Frequency of species', title='acutal species coverage')
plot2 <- ggplot(table, aes(x=position, y=allSpeciesCoverage, fill=species)) +
    geom_col() +
    theme(axis.text.x = element_text(angle=90, size=10)) +
    labs(x='position', y='Frequency of species', title='all species coverage')
filename <- tools::file_path_sans_ext(specieCoverage)
filename <- paste(filename, "allSpecies.png", sep=".")
ggsave(file=filename, arrangeGrob(plot1, plot2), width=24, height=12, units="cm", scale=2)

table_specieCoverage <- table %>% group_by(species) %>% summarise(speciesCount=sum(speciesCount), totalCount=sum(totalCount))
write.table(table_specieCoverage, file=paste(datadir, "speciesCoverage.tsv"), quote=FALSE, sep='\t', col.names = TRUE)

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

# library("tidyverse")
# library("gridExtra")
#
# ### load input files
# files <- list.files(path="/home/meijun/Documents/server/data/MSSPE-variant", recursive=TRUE, pattern="(?i)allOligo.set.coverage", full.name=TRUE)
#
# for (file in files){
#   table = read.delim(file=file, sep="\t", header=FALSE)
#   # names(table) = c("start", "end", "species", "speciesCount", "totalCount")
#   names(table) = c("start", "end", "chunk-start", "chunk-end", "species", "speciesCount", "totalCount")
#   table_specieCoverage <- table %>% group_by(species) %>% summarise(speciesCount=sum(speciesCount), totalCount=sum(totalCount))
#   table_specieCoverage$coverage <- table_specieCoverage$speciesCount/table_specieCoverage$totalCount
#   filename <- tools::file_path_sans_ext(file)
#   filename <- paste(filename, "speciesCoverage.png", sep=".")
#   plot <- ggplot(table_specieCoverage, aes(x=species, y=coverage)) +
#     geom_col() +
#     theme(axis.text.x = element_text(angle=90, size=10)) +
#     labs(x='species', y='Coverage of species', title='all species coverage')
#   ggsave(file=filename, plot, width=12, height=12, units="cm", scale=2)
# }
#
# table = read.delim(file=files[2], sep="\t", header=FALSE)
# names(table) = c("start", "end", "species", "speciesCount", "totalCount")
# table_specieCoverage <- table %>% group_by(species) %>% summarise(speciesCount=sum(speciesCount), totalCount=sum(totalCount))
# table_specieCoverage$coverage <- table_specieCoverage$speciesCount/table_specieCoverage$totalCount
# # filename <- tools::file_path_sans_ext(file)
# # filename <- paste(filename, "species.coverage.tsv", sep=".")
# # write.table(table_specieCoverage, file=filename, quote=FALSE, sep='\t', col.names = TRUE)
# plot <- ggplot(table_specieCoverage, aes(x=species, y=coverage)) +
#   geom_col() +
#   theme(axis.text.x = element_text(angle=90, size=10)) +
#   labs(x='species', y='Coverage of species', title='all species coverage')
# plot
