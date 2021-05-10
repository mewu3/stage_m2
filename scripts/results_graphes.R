library("tidyverse")
library("gridExtra")
library("reshape2")
library("seqinr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.6")

### load input files
datadir = "/home/meijun/Documents/server/data"
pipeline = c("MSSPE-basic", "MSSPE-variant")
virus = c("enterovirus", "coronavirus", "coronavirusWithoutSarsCov2", "enterovirusCurated", "coronavirusCurated", "coronavirusCuratedWithoutSarsCov2")
kmerSize = c("13", "15")
combinations = expand.grid(pipeline, virus, kmerSize)

for (row in 1:nrow(combinations)){
  pipeline = toString(combinations[row,1])
  virus = toString(combinations[row,2])
  kmerSize = toString(combinations[row,3])
  seqCount = 0 
  seqFile <- sprintf("%s/%s/%s/%s.uniq", datadir, pipeline, virus, virus)
  if ( (file.exists(seqFile)) & !str_detect(seqFile, "Curated") ) {
    commande <- sprintf("grep -c '^>' %s", seqFile)
    seqCount = system(commande, intern=TRUE)
  } 
  seqFile <- sprintf("%s/%s.fasta", datadir, virus)
  if ( (file.exists(seqFile)) & str_detect(seqFile, "Curated") ) {
    commande <- sprintf("grep -c '^>' %s", seqFile)
    seqCount = system(commande, intern=TRUE)
  }
  seqCount = as.numeric(seqCount)
  workdir <- sprintf("%s/%s/%s/kmer%s/evaluation", datadir, pipeline, virus, kmerSize)
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
    for (i in c("kmerCount.a0", "kmerCount.a1", "kmerCount.a2", "kmerCount.a3")){
      oligoCount[, i] <- oligoCount[, i]/seqCount * 100
      oligoCount[, i][oligoCount[, i]>100] <- 100
    }
    dfm <- pivot_longer(oligoCount, -position, names_to="states", values_to="value")
    plot <- ggplot(dfm, aes(x=position,y=value)) +
      geom_bar(aes(fill=states), stat = "identity", position = "dodge") +
      theme(axis.text.x = element_text(angle=90, size=9),
            legend.background=element_blank(),
            legend.title = element_blank()) +
      ylim(0,100) +
      xlab("Position") +
      ylab("Le comptage de kmer")
    filename <- paste(workdir, "kmerCount.png", sep="/")
    ggsave(plot, file=filename, width=30, height=12, units="cm", scale=2)
    
    ### species coverage ###
    finalSpeciesCoverage = sprintf("%s/allOligo.set.coverage", workdir)
    finalSpeciesCoverage = read.table(file=finalSpeciesCoverage,
                                      sep="\t",
                                      header=FALSE)
    names(finalSpeciesCoverage) = c("start", "end", "species", "speciesCount", "totalCount")
    finalSpeciesCoverage <- finalSpeciesCoverage %>% arrange(start, decreasing=TRUE)
    finalSpeciesCoverage$speciesCount[with(finalSpeciesCoverage, speciesCount > totalCount)] <- finalSpeciesCoverage$totalCount[with(finalSpeciesCoverage, speciesCount > totalCount)]
    finalSpeciesCoverage$position <- paste(finalSpeciesCoverage$start, finalSpeciesCoverage$end, sep="-")
    finalSpeciesCoverage$position <- factor(finalSpeciesCoverage$position, levels=unique(finalSpeciesCoverage$position))
    finalSpeciesCoverage[, "speciesCoverage"] = finalSpeciesCoverage[,"speciesCount"]/seqCount*100
    finalSpeciesCoverage[, "totalCoverage"] = finalSpeciesCoverage[,"totalCount"]/seqCount*100
    plot1 <- ggplot(finalSpeciesCoverage, aes(x=position, y=speciesCoverage, fill=species)) +
        geom_col() +
        theme(axis.text.x = element_text(angle=90, size=9),
              legend.background=element_blank(),
              legend.title = element_blank()) +
        ylim(0,100) +
        labs(x='Position', y='La fréquence de chaque espèce actuelle')
    totalCoverage = group_split(finalSpeciesCoverage, position)[[1]]
    plot2 <- ggplot(totalCoverage, aes(x=position, y=totalCoverage, fill=species)) +
        geom_col() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              legend.position="none",
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()) +
        labs(y="La fréquence de chaque espèce initiale") + 
        ylim(0,100) 
    # grid.arrange(plot2, plot1, widths=c(0.1,0.9), ncol=2)
    filename <- paste(workdir, "speciesCoverage.png", sep="/")
    ggsave(grid.arrange(plot2, plot1, widths=c(0.05,0.95), ncol=2), file=filename, width=30, height=12, units="cm", scale=2)
  }
}

# datadir <- sprintf("%s/%s/%s/kmer%s/evaluation", datadir, filePath, virus, kmerSize)
# 
# kmerCount <- list.files(path=datadir, pattern = "tsv", full.names = TRUE)
# specieCoverage = paste(datadir, "allOligo.set.coverage", sep="/")
# specieCoverage = paste(datadir, "allKmerCount.sorted.calculated.filtered.allOligo.set.coverage", sep="/")
# for (row in 1:nrow(combinations)){
#   pipeline = toString(combinations[row,1])
#   virus = toString(combinations[row,2])
#   kmerSize = toString(combinations[row,3])
#   seqCount = 0 
#   seqFile <- sprintf("%s/%s/%s/%s.uniq", datadir, pipeline, virus, virus)
#   if ( (file.exists(seqFile)) & !str_detect(seqFile, "Curated") ) {
#     commande <- sprintf("grep -c '^>' %s", seqFile)
#     seqCount = system(commande, intern=TRUE)
#   } 
#   seqFile <- sprintf("%s/%s.fasta", datadir, virus)
#   if ( (file.exists(seqFile)) & str_detect(seqFile, "Curated") ) {
#     commande <- sprintf("grep -c '^>' %s", seqFile)
#     seqCount = system(commande, intern=TRUE)
#   }
#   seqCount = as.numeric(seqCount)
#   workdir <- sprintf("%s/%s/%s/kmer%s/evaluation", datadir, pipeline, virus, kmerSize)
#   if (dir.exists(workdir)){
#     ### kmerCoverage counting ###
#     oligoCount_before <- sprintf("%s/allOligo_before.tsv", workdir)
#     oligoCount_before <- read.table(file=oligoCount_before,
#                                     sep="\t",
#                                     header=TRUE)
#     oligoCount_before <-oligoCount_before %>% select(start, end, kmerCount)
#     oligoCount_after1 <- sprintf("%s/allOligo_after1.tsv", workdir)
#     oligoCount_after1 <- read.table(file=oligoCount_after1,
#                                     sep="\t",
#                                     header=TRUE)
#     oligoCount_after1 <- oligoCount_after1 %>% select(start, end, kmerCount)
#     oligoCount_after2 <- sprintf("%s/allOligo_after2.tsv", workdir)
#     oligoCount_after2 <- read.table(file=oligoCount_after2,
#                                     sep="\t",
#                                     header=TRUE)
#     oligoCount_after2 <- oligoCount_after2 %>% select(start, end, kmerCount)
#     oligoCount_after3 <- sprintf("%s/allOligo_after3.tsv", workdir)
#     oligoCount_after3 <- read.table(file=oligoCount_after3,
#                                     sep="\t",
#                                     header=TRUE)
#     oligoCount_after3 <-oligoCount_after3 %>% select(start, end, kmerCount)
#     oligoCount <- left_join(x=oligoCount_before, y=oligoCount_after1, by=c("start", "end"), suffix=c(".a0", ".a1"))
#     oligoCount <- left_join(x=oligoCount, y=oligoCount_after2, by=c("start", "end"), suffix=c(".a1", ".a2"))
#     oligoCount <- left_join(x=oligoCount, y=oligoCount_after3, by=c("start", "end"), suffix=c(".a2", ".a3"))
#     oligoCount[is.na(oligoCount)] <- 0
#     oligoCount <- oligoCount %>% arrange(start, decreasing=TRUE)
#     oligoCount$position <- paste(oligoCount$start, oligoCount$end, sep="-")
#     oligoCount$position <- factor(oligoCount$position, levels=oligoCount$position)
#     oligoCount <- oligoCount[, !(names(oligoCount) %in% c("start", "end"))]
#     for (i in c("kmerCount.a0", "kmerCount.a1", "kmerCount.a2", "kmerCount.a3")){
#       oligoCount[, i] <- oligoCount[, i]/seqCount * 100
#       oligoCount[, i][oligoCount[, i]>100] <- 100
#     }
#     dfm <- pivot_longer(oligoCount, -position, names_to="variable", values_to="value")
#     plot <- ggplot(dfm, aes(x=position,y=value)) +
#       geom_bar(aes(fill=variable), stat = "identity", position = "dodge") +
#       theme(axis.text.x = element_text(angle=90, size=5)) +
#       ylim(0,100) +
#       xlab("position") +
#       ylab("kmer Coverage")
#     filename <- paste(workdir, "kmerCount.png", sep="/")
#     ggsave(plot, file=filename, width=30, height=12, units="cm", scale=2)
#     
#     ### species coverage ###
#     finalSpeciesCoverage = sprintf("%s/allOligo.set.coverage", workdir)
#     finalSpeciesCoverage = read.table(file=finalSpeciesCoverage,
#                                       sep="\t",
#                                       header=FALSE)
#     names(finalSpeciesCoverage) = c("start", "end", "species", "speciesCount", "totalCount")
#     finalSpeciesCoverage <- finalSpeciesCoverage %>% arrange(start, decreasing=TRUE)
#     finalSpeciesCoverage$speciesCount[with(finalSpeciesCoverage, speciesCount > totalCount)] <- finalSpeciesCoverage$totalCount[with(finalSpeciesCoverage, speciesCount > totalCount)]
#     finalSpeciesCoverage$position <- paste(finalSpeciesCoverage$start, finalSpeciesCoverage$end, sep="-")
#     finalSpeciesCoverage$position <- factor(finalSpeciesCoverage$position, levels=unique(finalSpeciesCoverage$position))
#     finalSpeciesCoverage[, "speciesCoverage"] = finalSpeciesCoverage[,"speciesCount"]/seqCount*100
#     finalSpeciesCoverage[, "totalCoverage"] = finalSpeciesCoverage[,"totalCount"]/seqCount*100
#     plot1 <- ggplot(finalSpeciesCoverage, aes(x=position, y=speciesCoverage, fill=species)) +
#       geom_col() +
#       theme(axis.text.x = element_text(angle=90, size=5)) +
#       ylim(0,100) +
#       labs(x='position', y='Frequency of species', title='acutal species coverage')
#     totalCoverage = group_split(finalSpeciesCoverage, position)[[1]]
#     plot2 <- ggplot(totalCoverage, aes(x=position, y=totalCoverage, fill=species)) +
#       geom_col() +
#       theme(axis.title.x=element_blank(),
#             axis.text.x=element_blank(),
#             axis.ticks.x=element_blank(),
#             legend.position="none",
#             panel.border = element_blank(),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(),
#             panel.background = element_blank()) +
#       ylim(0,100) 
#     # grid.arrange(plot2, plot1, widths=c(0.1,0.9), ncol=2)
#     filename <- paste(workdir, "speciesCoverage.png", sep="/")
#     ggsave(grid.arrange(plot2, plot1, widths=c(0.1,0.9), ncol=2), file=filename, width=30, height=12, units="cm", scale=2)
#     
#   }
  
### global variables
# totolSeq = 3265 # enterovirus
# totolSeq = 75959 # coronavirus
# totolSeq = 119 # enterovirusCurated
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
# table = read_delim(file=specieCoverage, delim="\t", col_names=FALSE)
# names(table) = c("start", "end", "chunk-start", "chunk-end", "species", "speciesCount", "totalCount")
# table[, "speciesCoverage"] = table[,"speciesCount"]/totolSeq*100
# table[, "allSpeciesCoverage"] = table[,"totalCount"]/totolSeq*100
# table <- arrange(table, start)
# table[, "position"] <- paste(table$`chunk-start`, table$`chunk-end`, sep="-")
# table$position <- factor(table$position, levels=unique(table$position))

# all together
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
# 
# table_specieCoverage <- table %>% group_by(species) %>% summarise(speciesCount=sum(speciesCount), totalCount=sum(totalCount))
# write.table(table_specieCoverage, file=paste(datadir, "speciesCoverage.tsv"), quote=FALSE, sep='\t', col.names = TRUE)

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
