library("tidyverse")
library("gridExtra")
library("reshape2")
library("plotrix")
library("svglite")

workdir <- "/home/meijun/Documents/server/data/MSSPE-variant/enterovirusCurated/kmer13/evaluation"
seqFile <- "/home/meijun/Documents/server/data/MSSPE-variant/enterovirusCurated/enterovirusCurated.msa"
commande <- sprintf("grep -c '^>' %s", seqFile)
seqCount <- system(commande, intern = TRUE)
seqCount <- as.numeric(seqCount)

genomeCoverage = sprintf("%s/genomeCoverage.tsv", workdir)
genomeCoverage = read.table(file=genomeCoverage, sep="\t", header=FALSE)
names(genomeCoverage) = c("Oligonucleotides", "Espece", "Espece_reel", "Espece_total")

genomeCoverage$Oligonucleotides <- factor(genomeCoverage$Oligonucleotides, levels = unique(genomeCoverage$Oligonucleotides))

genomeCoverage$Couverture = genomeCoverage$Espece_reel / seqCount * 100

plot <- ggplot(genomeCoverage, aes(x=Oligonucleotides, y=Couverture, fill=Espece)) +
                geom_col() +
                labs(title="Pourcentage de génomes ciblés par au moins x oligonucléotides", fill="Espèces", x="Nombre d'oligonucléotides", y="Pourcentage de génomes") +
                theme(plot.title = element_text(size=15, face="bold"),
                      axis.title = element_text(size = 15),
                      axis.text.x = element_text(size = 10),
                      axis.text.y = element_text(size = 10),
                      legend.title = element_text(size=15),
                      legend.text = element_text(size=13),
                      legend.background = element_blank()) +
                # scale_x_continuous("Oligonucleotides", labels=as.character(Oligonucleotides), breaks=Oligonucleotides) + 
                ylim(0, 100) 

filename <- paste(workdir, "genomeCoverage.png", sep = "/")
ggsave(plot, file = filename, width = 30, height = 15, units = "cm")

