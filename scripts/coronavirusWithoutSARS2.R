library("tidyverse")
library("gridExtra")
library("reshape2")

### load input files
datadir = "/datater/wu/data/pre-historic"
pipeline = c("MSSPE-basic", "MSSPE-variant")
virus = c("coronavirusWithoutSarsCov2", "coronavirusCuratedWithoutSarsCov2")
kmerSize = c("13", "15")
combinations = expand.grid(pipeline, virus, kmerSize)

for (row in 1:nrow(combinations)){
  pipeline = toString(combinations[row, 1])  
  virus = toString(combinations[row,2])
  kmerSize = toString(combinations[row,3])
  file <- sprintf("%s/%s_%s_kmer%s_statistic.out", datadir, pipeline, virus, kmerSize)
  if (file.exists(file)){
      table <- read.table(file=file, sep="\t", header=TRUE)
      table[, "stricOligoCount"] <- table[, "counts_x"]/table[, "SeqCount"]
      table[, "notStricOligoCount"] <- table[, "counts_y"]/table[, "SeqCount"]
      table <- table %>% arrange(start, decreasing=TRUE)
      table$position <- paste(table$start, table$end, sep="-")
      table$position <- factor(table$position, levels=table$position)
      table <- table[, !(names(table) %in% c("start", "end", "counts_x", "counts_y", "SeqCount"))]
      table[, "stricOligoCount"] = table[,"stricOligoCount"] * 100
      table[,"notStricOligoCount"] = table[,"notStricOligoCount"]*100
      
      print(file)
      print(mean(table[,"notStricOligoCount"]))
#       dfm <- pivot_longer(table, -position, names_to="variable", values_to="value")
#       plot <- ggplot(dfm, aes(x=position,y=value)) +
#         geom_bar(aes(fill=variable), stat = "identity", position = "dodge") + 
#         theme(axis.text.x = element_text(angle=90, size=5)) +
#         xlab("position") + 
#         ylab("%Seq cible") + 
#         ylim(0,100)
#       filename <- tools::file_path_sans_ext(file)
#       filename <- paste(filename, ".png", sep="")
#       ggsave(plot, file=filename)
  }    
}

