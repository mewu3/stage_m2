library("tidyverse")
library("gridExtra")
library("reshape2")
library("seqinr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.6")
library("gg.gap", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.6")
library("plotrix")
library("svglite")

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
  
  if (str_detect(pipeline, "basic")) {
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
  }
  
  if (str_detect(pipeline, "variant")) {
    seqFile <- sprintf("%s/%s/%s/%s0.99.msa", datadir, pipeline, virus, virus)
    if ( (file.exists(seqFile)) & str_detect(seqFile, "variant") ) {
      commande <- sprintf("grep -c '^>' %s", seqFile)
      seqCount = system(commande, intern=TRUE)
    }
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
    oligoCount <- oligoCount %>% arrange(end, decreasing=TRUE)
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
    # print(plot)

    ### species coverage per position ###
    finalSpeciesCoverage = sprintf("%s/allOligo.set.coverage", workdir)
    finalSpeciesCoverage = read.table(file=finalSpeciesCoverage,
                                      sep="\t",
                                      header=FALSE)
    names(finalSpeciesCoverage) = c("start", "end", "species", "speciesCount", "totalCount")

    finalSpeciesCoverage <- finalSpeciesCoverage %>% filter(totalCount !=0)

    finalSpeciesCoverage$coverage <- finalSpeciesCoverage$speciesCount/finalSpeciesCoverage$totalCount
    finalSpeciesCoverage$species <- str_replace_all(finalSpeciesCoverage$species, "Betacoronavirus 1", "OC43")
    finalSpeciesCoverage$species <- str_replace_all(finalSpeciesCoverage$species, "Human coronavirus HKU1", "HKU1")
    finalSpeciesCoverage$species <- str_replace_all(finalSpeciesCoverage$species, "Human coronavirus 229E", "229E")
    finalSpeciesCoverage$species <- str_replace_all(finalSpeciesCoverage$species, "Human coronavirus NL63", "NL63")
    finalSpeciesCoverage$species <- str_replace_all(finalSpeciesCoverage$species, "Middle East respiratory syndrome-related coronavirus", "MERS")
    finalSpeciesCoverage$species <- str_replace_all(finalSpeciesCoverage$species, "Severe acute respiratory syndrome-related coronavirus", "SARS")

    finalSpeciesCoverage <- finalSpeciesCoverage %>% arrange(end, decreasing=TRUE)
    finalSpeciesCoverage$speciesCount[with(finalSpeciesCoverage, speciesCount > totalCount)] <- finalSpeciesCoverage$totalCount[with(finalSpeciesCoverage, speciesCount > totalCount)]
    finalSpeciesCoverage$position <- paste(finalSpeciesCoverage$start, finalSpeciesCoverage$end, sep="-")
    finalSpeciesCoverage$position <- factor(finalSpeciesCoverage$position, levels=unique(finalSpeciesCoverage$position))

    seqFile <- sprintf("%s/%s/%s0.99", datadir, pipeline, virus)
    if ( (file.exists(seqFile)) & str_detect(pipeline, "variant") ) {
      commande <- sprintf("grep -c '^>' %s", seqFile)
      seqCount = system(commande, intern=TRUE)
    }

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
    filename <- paste(workdir, "speciesCoverage.png", sep="/")
    ggsave(grid.arrange(plot2, plot1, widths=c(0.05,0.95), ncol=2), file=filename, width=30, height=12, units="cm", scale=2)
    # print(plot1)

    ### species coverage in total ###
    totalSpeciesCoverage <- finalSpeciesCoverage %>%
      group_by(species) %>%
      summarise(speciesCount=sum(speciesCount), totalCount=sum(totalCount))
    totalSpeciesCoverage$coverage <- totalSpeciesCoverage$speciesCount/totalSpeciesCoverage$totalCount *100
    plot3 <- ggplot(totalSpeciesCoverage, aes(x=species, y=coverage, fill=species)) +
      geom_col() +
      theme(axis.text.x = element_text(angle=90, size=9)) +
      ylim(0,100)
    # plot3 <- gg.gap(
    #   plot=plot3,
    #   segments=c(60,90),
    #   tick_width = 10,
    #   rel_heights = c(0.9, 0, 0.1),
    #   ylim=c(0,100)
    # )
    filename <- paste(workdir, "totalSpeciesCoverage.png", sep="/")
    ggsave(plot3, file=filename, scale=2)
    # print(plot3)
  }
}
