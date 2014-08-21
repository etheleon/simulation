#!/usr/bin/env Rscript

library(Metamaps2)

data=read.table("../data/topChosenGenera.txt", sep="\t",h=T)

query="MATCH (genus:genus)
WHERE genus.name = {taxaname} 
WITH genus 
MATCH 
 path=genus-[:childof*]->(king:superkingdom)
RETURN
 genus.name as NAME, 
 genus.taxid as TAXID,
 last(extract(n in nodes(path)| n.name)) AS SUPERKINGDOM"

taxondf=do.call(rbind,lapply(unique(data$taxon), function(x) {
	dbquery(query=query, param=list(taxaname=x), cypherurl = "metamaps.scelse.nus.edu.sg:7474/db/data/cypher")
}))
taxondf =setNames(as.data.frame(t(apply(taxondf,1,unlist))),c("name","taxid","superkingdom"))	#problematic df structure

taxondf=subset(taxondf, superkingdom %in% c("Bacteria", "Archaea"))

write.table(taxondf, file="out/sim.0020.out.txt",quote=F,sep="\t",row.names=F)
