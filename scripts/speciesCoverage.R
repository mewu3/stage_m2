library("tidyverse")
library("gridExtra")
library("reshape2")
library("seqinr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.6")
library("gg.gap", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.6")
library("plotrix")
library("svglite")


finalSpeciesCoverage = read.table(file=finalSpeciesCoverage,
                                  sep="\t",
                                  header=FALSE)
names(finalSpeciesCoverage) = c("start", "end", "species", "speciesCount", "totalCount")

totalSpeciesCoverage$coverage <- totalSpeciesCoverage$speciesCount/totalSpeciesCoverage$totalCount *100
plot3 <- ggplot(totalSpeciesCoverage, aes(x=species, y=coverage, fill=species)) +
  geom_col() +
  theme(axis.text.x = element_text(angle=90, size=9)) +
  ylim(0,100)
