#!/usr/bin/env Rscript

library(tidyverse)

file = "allOligo_before.tsv"
# file = "allOligo_after1.tsv"
# file = "allOligo_after2.tsv"
# file = "allOligo_after3.tsv"
filePath = sprintf("/home/meijun/Documents/result/%s", file)

table = read_tsv(file = filePath, 
                 col_names = TRUE)

ggplot(table, aes(position, kmerCount)) +
  geom_col() +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(breaks=table$position)
  ylim(0, 3265)
