library(tidyverse)

datadir <- "/home/meijun/Documents/server/data/enterovirus/evaluation13"

files <- list.files(path=datadir, pattern='.tsv', full.names = TRUE)
specieCoverage = paste(datadir, "allOligos_reverse.set.coverage", sep="/")

totolSeq = 3265

savePlot <- function(myPlot, plotname) {
  pdf(plotname)
  print(myPlot)
  dev.off()
}

### the most abondunce kmer coverage 
for (file in files) {
  filename = tools::file_path_sans_ext(file)
  table = read.delim(file=file, sep="\t")
  table[,2] <- table[,2]/totolSeq *100
  table2 <- table %>% separate(position, into = c("position1", "position2"))
  table3 <- cbind(table$position, table2)
  cols.num <- c("position1", "position2")
  table3[cols.num] <- sapply(table3[cols.num],as.numeric)
  colnames(table3)[1] <- "position"
  table3 <- table3 %>%
    arrange(position1)
  table3$position <- factor(table3$position, levels=table3$position)
  plot <- ggplot(table3, aes(x=position, y=kmerCount)) +
    geom_col() +
    theme(axis.text.x = element_text(angle=90, size=10)) +
    ylim(0,100) +
    ylab("kmer Coverage")
  plotname = paste(filename, ".pdf", sep="")
  savePlot(plot, plotname)
}

### species coverage 
table = read.delim(file=specieCoverage, sep="\t", header=FALSE)
names(table) = c("position", "species", "count", "totalCount")
table[,3] <- table[,3]/table[,4] *100
table2 <- table %>% separate(position, into = c("position1", "position2"))
table3 <- cbind(table$position, table2)
cols.num <- c("position1", "position2")
table3[cols.num] <- sapply(table3[cols.num],as.numeric)
colnames(table3)[1] <- "position"
table3 <- table3 %>%
  arrange(position1)
table3$position <- factor(table3$position, levels=unique(table3$position))

lsspecies = unique(table3$species)

for (spec in lsspecies){
  speTable = table3 %>% filter(species == spec)
  plot <- ggplot(speTable, aes(x=position, y=count, fill=species)) +
    geom_col() +
    theme(axis.text.x = element_text(angle=90, size=10)) + 
    ylim(0,100) + 
    labs(x='position', y='Frequency of species')
  filename = paste(datadir, spec, sep="/")
  plotname = paste(filename, ".pdf", sep="")
  savePlot(plot, plotname)
}

ggplot(table3, aes(x=position, y=count, fill=species)) +
  geom_col() +
  theme(axis.text.x = element_text(angle=90, size=10)) + 
  # ylim(0,100) + 
  labs(x='position', y='Frequency of species')
plotname = paste(datadir, "allSpecies.pdf", sep="/")
savePlot(plot, plotname)
