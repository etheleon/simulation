#!/usr/bin/env Rscript

library(Metamaps2)
library(dplyr)
#Using prokaryotes
#Description
#prokaryotes.txt: 	Prokaryotic genome sequencing projects
#			excluding projects that represent only plasmids
data=read.table("data/genomereports/prokaryotes.txt", sep="\t",fill=T,comment.char="",h=T)

#Choose only genera with refseq sequences
refseq.data = subset(data, Chromosomes.RefSeq != "-")

#get the genus of these sequences
query="START basetaxa=node:ncbitaxid(taxid={taxid})
MATCH basetaxa-[:childof*0..]->(higher)
WHERE higher:genus OR higher:family OR higher:order OR higher:class OR higher:phylum
RETURN head(labels(higher)), higher.taxid,basetaxa.taxid"

#match genera
taxondf=	setNames(do.call(rbind,lapply(unique(refseq.data$TaxID), function(x){
	dbquery(
		query=query,
	    	cypherurl="metamaps.scelse.nus.edu.sg:7474/db/data/cypher", 
	    	params=list(taxid=x)	
	    	)
})), c("Rank", "taxid", "origin"))
taxondf =as.data.frame(t(apply(taxondf,1,unlist)))	#problematic df structure
taxondf.genus=taxondf %.% group_by(origin) %.% filter(Rank == 'genus')

#the Taxid of the observable taxa in the sludge
abu=read.table("data/top500genera.gDNA",sep="\t",h=T)

query1=	"MATCH (genus:genus) 
WHERE genus.name = {nameoftaxa}
return genus.name, genus.taxid"

abu2=setNames(do.call(rbind,lapply(as.character(abu$taxon),function(x){
    dbquery(
    	query=query1, 
    	params=list(nameoftaxa = x), 
    	cypherurl="metamaps.scelse.nus.edu.sg:7474/db/data/cypher"
	)
    }))
, c("taxon","taxid"))

abu=merge(abu,abu2, by="taxon")
abu=abu[order(abu$rank),]

#Combined genome annotation df with the genera annotation
combined=merge(taxondf.genus, abu, by.x="taxid", by.y="taxid")
combined=combined[order(combined$rank),]

#Finds how 1:nth genera to consider 
u=100;i=1
while(i<100){
    	i <- sum(1:u %in% combined$rank)
	u <- u + 1
}

#sum(1:u %in% subset(combined, rank <= u)$rank)	#100 unique genera

combined2=subset(combined,rank<=u)
combined3=merge(combined2, refseq.data, by.x="origin", by.y="TaxID")

#Choosen genome::luckydraw
set.seed(5)
chosen_genomes=do.call(rbind,lapply(unique(combined3$taxid), function(x) { 
genus=subset(combined3, taxid == x)
genus[sample(1:nrow(genus),1),]
}))

#OUTPUT

#1 for parsing the genomes
write.table(chosen_genomes[,c("origin","taxid","taxon","Chromosomes.RefSeq")], file="out/sim.0010.out.txt", row.names=F, sep="\t", quote=F)

#2 for removing the related sequence
chosenones=unique(combined3[,c("taxon","total","rank")])
chosenones=chosenones[order(chosenones$rank),]
write.table(chosenones, file="data/topChosenGenera.txt", quote=F, row.names=F,sep="\t")
