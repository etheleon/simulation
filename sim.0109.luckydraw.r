#!/usr/bin/env Rscript

#Aim: chooses one gi with refseq sequence for each taxid
data=setNames(read.table("out/sim.0108.out.txt",sep="\t",h=F), c("gi","base","genus"))

new.data=
do.call(rbind,
lapply(unique(data$genus), function(x) { 
    data2 = subset(data, genus == x)
    data2[sample(nrow(data2),1),]
    })
)

#-- only 94 have a refseq sequence
write.table(new.data, file="out/sim.0109.out.txt",quote=F,row.names=F)

library(RCurl)
library(RJSONIO)
idk=dbquery(
#query = "START taxon=node:ncbitaxid(taxid={taxa}) RETURN taxon.taxid,taxon.name",
query = "START inputko=node:koid(ko={koid}) MATCH inputko--(cpd:`cpd`)--(outputko:`ko`) return inputko.ko, outputko.ko",
	params = list(koid= "ko:K00020"), 
	cypherurl = "192.168.100.1:7474/db/data/cypher"
	)
