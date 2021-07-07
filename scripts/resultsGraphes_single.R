library("tidyverse")
library("gridExtra")
library("reshape2")
library("seqinr")
library("plotrix")
library("svglite", lib = "/home/meijun/Programmes/miniconda3/lib/R/library/")
# library("seqinr", lib="/datater/wu/miniconda3/lib/R/library/")
# library("plotrix", lib="/datater/wu/miniconda3/lib/R/library/")
# library("svglite", lib="/datater/wu/miniconda3/lib/R/library/")

### load input files ###
workdir <- "/home/meijun/Documents/server/data/test_clustering/enterovirus_normalized/kmer13/evaluation"
seqFile <- "/home/meijun/Documents/server/data/test_clustering/enterovirus_normalized.fasta"
commande <- sprintf("grep -c '^>' %s", seqFile)
seqCount <- system(commande, intern = TRUE)
seqCount <- as.numeric(seqCount)

### oligoCount ###
oligoCount_before <- sprintf("%s/allOligo_before.tsv", workdir)
oligoCount_before <- read.table(file = oligoCount_before,
                                sep = "\t",
                                header = TRUE)
oligoCount_before <-oligoCount_before %>% select(start, end, kmerCount)

oligoCount_after1 <- sprintf("%s/allOligo_after1.tsv", workdir)
oligoCount_after1 <- read.table(file = oligoCount_after1,
                                sep = "\t",
                                header = TRUE)
oligoCount_after1 <- oligoCount_after1 %>% select(start, end, kmerCount)

oligoCount_after2 <- sprintf("%s/allOligo_after2.tsv", workdir)
oligoCount_after2 <- read.table(file = oligoCount_after2,
                                sep = "\t",
                                header = TRUE)
oligoCount_after2 <-oligoCount_after2 %>% select(start, end, kmerCount)

oligoCount_after3 <-sprintf("%s/allOligo_after3.tsv", workdir)
oligoCount_after3 <- read.table(file = oligoCount_after3,
                                sep = "\t",
                                header = TRUE)
oligoCount_after3 <-oligoCount_after3 %>% select(start, end, kmerCount)

oligoCount <- left_join(x = oligoCount_before,
                        y = oligoCount_after1,
                        by = c("start", "end"),
                        suffix = c(".a0", ".a1"))
oligoCount <- left_join(x = oligoCount,
                        y = oligoCount_after2,
                        by = c("start", "end"),
                        suffix = c(".a1", ".a2"))
oligoCount <- left_join(x = oligoCount,
                        y = oligoCount_after3,
                        by = c("start", "end"),
                        suffix = c(".a2", ".a3"))

oligoCount[is.na(oligoCount)] <- 0

oligoCount <- oligoCount %>% arrange(end, decreasing = TRUE)
oligoCount$position <- paste(oligoCount$start, oligoCount$end, sep = "-")
oligoCount$position <- factor(oligoCount$position, levels = oligoCount$position)
oligoCount <- oligoCount[,!(names(oligoCount) %in% c("start", "end"))]

for (i in c("kmerCount.a0", "kmerCount.a1", "kmerCount.a2", "kmerCount.a3")) {
    oligoCount[, i] <- oligoCount[, i] / seqCount * 100
    oligoCount[, i][oligoCount[, i] > 100] <- 100
}

names(oligoCount)[names(oligoCount) == "kmerCount.a0"] <- "Sélection 1"
names(oligoCount)[names(oligoCount) == "kmerCount.a1"] <- "Filtrage 1"
names(oligoCount)[names(oligoCount) == "kmerCount.a2"] <- "Filtrage 2"
names(oligoCount)[names(oligoCount) == "kmerCount.a3"] <- "Sélection 2"

dfm <- pivot_longer(oligoCount, -position, names_to = "states", values_to = "value")

dfm$states <- factor(dfm$states, levels = unique(dfm$states))

plot <- ggplot(dfm, aes(x = position, y = value)) +
        geom_bar(aes(fill = states),
                 stat = "identity",
                 position = "dodge") +
        theme(axis.text.x = element_text(angle = 90, size = 10),
              axis.text.y = element_text(size = 10),
              axis.title = element_text(size = 15, face = "bold"),
              legend.background = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(size=13)) +
        # scale_x_discrete(breaks = levels(dfm$position)[c(T, rep(F,5))] ) +
        ylim(0, 100) +
        xlab("Position") +
        ylab("Couverture des k-mers (%)")
filename <- paste(workdir, "kmerCount.png", sep = "/")
ggsave(plot, file = filename, width = 30, height = 15, units = "cm")

plot <- ggplot(oligoCount, aes(x =position, y =`Sélection 2`)) +
        geom_bar(stat="identity", aes(fill="Sélection 2")) +
        theme(axis.text.x = element_text(angle = 90, size = 10),
              axis.text.y = element_text(size = 10),
              axis.title = element_text(size = 15, face = "bold"),
              legend.background = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(size=13)) +
        # scale_x_discrete(breaks = levels(dfm$position)[c(T, rep(F,5))] ) +
        scale_fill_manual(breaks = 'Sélection 2', values = '#C77CFF') +
        ylim(0, 100) +
        xlab("Position") +
        ylab("Couverture des k-mers (%)")
filename <- paste(workdir, "kmerCountFinal.png", sep = "/")
ggsave(plot, file = filename, width = 30, height = 15, units = "cm")

# if (str_detect(pipeline, "variant")) {
#     seqFile <-
#         sprintf("%s/%s/%s/%s.msa", datadir, pipeline, virus, virus)
#     commande <- sprintf("grep -c '^>' %s", seqFile)
#     seqCount = system(commande, intern = TRUE)
# }
# seqCount = as.numeric(seqCount)

### species coverage per position ###
finalSpeciesCoverage = sprintf("%s/allOligo.set.coverage", workdir)
finalSpeciesCoverage = read.table(file = finalSpeciesCoverage,
                                  sep = "\t",
                                  header = FALSE)
names(finalSpeciesCoverage) = c("start", "end", "species", "speciesCount", "totalCount")

finalSpeciesCoverage <- finalSpeciesCoverage %>% filter(totalCount != 0)

finalSpeciesCoverage$coverage <- finalSpeciesCoverage$speciesCount / finalSpeciesCoverage$totalCount
finalSpeciesCoverage$species <- str_replace_all(finalSpeciesCoverage$species, "Betacoronavirus 1", "OC43")
finalSpeciesCoverage$species <- str_replace_all(finalSpeciesCoverage$species, "Human coronavirus HKU1", "HKU1")
finalSpeciesCoverage$species <- str_replace_all(finalSpeciesCoverage$species, "Human coronavirus 229E", "229E")
finalSpeciesCoverage$species <- str_replace_all(finalSpeciesCoverage$species, "Human coronavirus NL63", "NL63")
finalSpeciesCoverage$species <- str_replace_all(finalSpeciesCoverage$species, "Middle East respiratory syndrome-related coronavirus", "MERS" )
finalSpeciesCoverage$species <- str_replace_all(finalSpeciesCoverage$species, "Severe acute respiratory syndrome-related coronavirus", "SARS")

finalSpeciesCoverage <- finalSpeciesCoverage %>% arrange(end, decreasing = TRUE)
finalSpeciesCoverage$speciesCount[with(finalSpeciesCoverage, speciesCount > totalCount)] <- finalSpeciesCoverage$totalCount[with(finalSpeciesCoverage, speciesCount > totalCount)]
finalSpeciesCoverage$position <- paste(finalSpeciesCoverage$start, finalSpeciesCoverage$end, sep = "-")
finalSpeciesCoverage$position <- factor(finalSpeciesCoverage$position, levels = unique(finalSpeciesCoverage$position))

finalSpeciesCoverage[, "speciesCoverage"] = finalSpeciesCoverage[, "speciesCount"] /seqCount * 100
finalSpeciesCoverage[, "totalCoverage"] = finalSpeciesCoverage[, "totalCount"] / seqCount * 100

plot1 <- ggplot(finalSpeciesCoverage, aes(x = position, y = speciesCoverage, fill = species)) +
         geom_col() +
         labs(title="Couverture par les oligonucléotides", fill="Espèces") + 
         theme(axis.text.x = element_text(angle = 90, size = 10),
               axis.text.y = element_text(size = 10),
               axis.title = element_text(size = 15, face="bold"),
               legend.background = element_blank(),
               plot.title = element_text(size=15, face="bold"),
               legend.title = element_text(size=15, face="bold"),
               legend.text = element_text(size=13)) +
        # scale_x_discrete(breaks = levels(dfm$position)[c(T, rep(F,5))] ) +
        ylim(0, 100) +
        labs(x = 'Position', y = "Couverture d'espèces (%)")

totalCoverage = group_split(finalSpeciesCoverage, position)[[1]]
plot2 <-ggplot(totalCoverage, aes(x = position, y = totalCoverage, fill = species)) +
        geom_col() +
        labs(y = "Abondance en espèces (%)") + 
        theme(axis.title.x = element_blank(),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             legend.position = "none",
             axis.title = element_text(size = 15, face = "bold"),
             panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank()) +
        ylim(0, 100)

filename <- paste(workdir, "speciesCoverage.png", sep = "/")
ggsave(grid.arrange(plot2, plot1, widths = c(0.1, 0.9), ncol = 2), file = filename, width = 30, height = 15, units = "cm")

### species coverage in total ###
totalSpeciesCoverage <- finalSpeciesCoverage %>% group_by(species) %>% summarise(speciesCount = sum(speciesCount), totalCount = sum(totalCount))
totalSpeciesCoverage$coverage <- totalSpeciesCoverage$speciesCount / totalSpeciesCoverage$totalCount * 100

plot3 <- ggplot(totalSpeciesCoverage, aes(x = species, y = coverage, fill = species)) +
         geom_col() +
         labs(title="Couverture par espèces", fill="Espèces", y="Couverture (%)") +
         theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.text.y = element_text(size = 13),
               plot.title = element_text(size=15, face="bold"),
               axis.title = element_text(size = 15, face = "bold"),
               legend.title = element_text(size=15, face="bold"),
               legend.text = element_text(size=13)) +
         geom_text(aes(label=round(coverage, digits =1)), vjust=-0.5) + 
         ylim(0, 100)

filename <-paste(workdir, "totalSpeciesCoverage.png", sep = "/")
ggsave(plot3, file = filename, width = 15, height = 15, units = "cm")

