#!/usr/bin/env Rscript

library(RCurl)
library(RJSONIO)
library(dplyr)
cypherurl = "192.168.100.1:7474/db/data/cypher"

#AIM: Just take 500 
dbquery<-function(query, cypherurl, params){
df_raw = fromJSON(
getURL(
	cypherurl, 
    	customrequest='POST', 
        httpheader=c('Content-Type'='application/json'),
        postfields=toJSON(list(query=query, params=params))
              )	
)
as.data.frame(do.call(rbind,lapply(df_raw$data, function(x) matrix(x, nrow=1)))) 
}

#Using prokaryotes
#Description
#prokaryotes.txt: 	Prokaryotic genome sequencing projects
#			excluding projects that represent only plasmids
data=read.table("data/genomereports/prokaryotes.txt", sep="\t",fill=T,comment.char="",h=T)

#Choose only genera with refseq sequences
refseq.data = filter(data, Chromosomes.RefSeq != "-", Status == 'Gapless Chromosome')

#get the genus of these sequences
query="START basetaxa=node:ncbitaxid(taxid={taxid})
MATCH basetaxa-[:childof*0..]->(higher)
WHERE higher:genus OR higher:family OR higher:order OR higher:class OR higher:phylum
RETURN head(labels(higher)), higher.taxid,basetaxa.taxid"

#match genera
taxondf=	setNames(do.call(rbind,lapply(unique(refseq.data$TaxID), function(x){
dbquery(
		query=query,
	    	cypherurl=cypherurl, 
	    	params=list(taxid=x)	
	    	)
})), c("Rank", "taxid", "origin"))
taxondf =as.data.frame(t(apply(taxondf,1,unlist)))	#problematic df structure
taxondf.genus=taxondf %.% group_by(origin) %.% filter(Rank == 'genus')

#the Taxid of the observable taxa in the sludge
abu=read.table("data/top500genera.gDNA",sep="\t",h=T)

#Might have more than 1 taxa with the same name, 
query1="MATCH (genus:genus)
WHERE genus.name = {taxaname} 
WITH genus 
MATCH 
 path=genus-[:childof*]->(king:superkingdom)
RETURN
 genus.name as NAME, 
 genus.taxid as TAXID,
 last(extract(n in nodes(path)| n.name)) AS SUPERKINGDOM"

abu2=setNames(do.call(rbind,lapply(as.character(abu$taxon),function(x){
    dbquery(
    	query=query1, 
    	params=list(taxaname= x), 
    	cypherurl=cypherurl
	)
    }))
, c("taxon","taxid","superkingdom"))
abu2=subset(abu2, superkingdom %in% c("Bacteria", "Archaea"))
#genera from bacteria/Archaea in top500 genera 
#413

abu2=abu2[,c(1,2)]

#for writing to output as input to sim.0200
abuu=setNames(as.data.frame(do.call(cbind,
lapply(1:ncol(abu2), function(x) {
unlist(abu2[,x])
    }
    ))), colnames(abu2)
)
abuu = merge(abu, abuu, by="taxon")
write.table(abuu, file="out/sim.0101.out2.txt",quote=F, row.names=F,sep="\t")


abu=merge(abu,abu2, by="taxon")
abu=abu[order(abu$rank),]

#Combined genome annotation df with the genera annotation
combined=merge(taxondf.genus, abu, by.x="taxid", by.y="taxid")
combined=combined[order(combined$rank),]

#Finds how 1:nth genera to consider 
#dun need right
#u=100;i=1
#while(i<100){
#    	i <- sum(1:u %in% combined$rank)
#	u <- u + 1
#}

#sum(1:u %in% subset(combined, rank <= u)$rank)	#100 unique genera
#combined2=subset(combined,rank<=u)

#264 genera with complete genomes
combined3=merge(combined, refseq.data, by.x="origin", by.y="TaxID")

#Choosen genome::luckydraw
set.seed(5)
chosen_genomes=do.call(rbind,lapply(unique(combined3$taxid), function(x) { 
genus=subset(combined3, taxid == x)
genus[sample(1:nrow(genus),1),]
}))

#OUTPUT
#for genera with complete genomes
write.table(chosen_genomes[,c("origin","taxid","taxon","Chromosomes.RefSeq")], 
file="out/sim.0101.out.txt", row.names=F, sep="\t", quote=F)

#2 for removing the related sequence
chosenones=unique(combined3[,c("taxon","total","rank")])
chosenones=chosenones[order(chosenones$rank),]
write.table(chosenones, file="data/topChosenGenera.txt", quote=F, row.names=F,sep="\t")

#3 To process those without
write.table(
x=data.frame(
taxid=
do.call(c,abu2$taxid[!abu2$taxid %in% unique(combined$taxid)])), 
sep="\t", 
row.names=F, 
quote=F,
file="out/sim.0101.missing.txt")
