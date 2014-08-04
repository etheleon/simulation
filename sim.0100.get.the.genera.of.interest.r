#!/usr/bin/env Rscript

library(dplyr)
load("data/Abundance_Function/mRNA.1000.genus.rda")

genera=genus %.%
filter(taxon != 'Unclassified') %.%
group_by(taxon) %.% 
summarise(total=sum(c1.mean)) %.% 
arrange(desc(total))

top100genera = cbind(genera[1:100,], data.frame(rank=1:100))
write.table(top100genera, file="data/top100genera.gDNA",sep="\t",row.names=F,quote=F)
