#!/usr/bin/env Rscript

library(RCurl)
library(data.table)
library(RJSONIO)
library(dplyr)
library(ggplot2)

#Reads in the output from the initial parent-child link and the gi-child file
parentchild=read.table("out/sim.0101.out.txt", h=T, sep="\t", skip=2)
removedgenomes=read.table("out/sim.0107.out.txt",h=T,sep="\t")
combined = merge(parentchild, removedgenomes, by.x="targetID", by.y="taxid")

load("~/abu_genes/data/Abundance_Function/g454.1000.genus.rda")
top500=head(tail(genus.454 %.% group_by(taxon) %.% summarise(abundance=sum(g1.mean)) %.% arrange(desc(abundance)) , n=-1), n=500)

#
summaryTaxa = combined %.%
group_by(originID,targetID) %.%
summarise(sizePerChild = n())

#
summaryTaxa0 = combined %.% 
group_by(originID) %.%
summarise(numUniqueChild = length(unique(targetID)))

summaryTaxa1=merge(summaryTaxa, summaryTaxa0, by="originID")

source("/export2/home/uesu/github/altimit/metamaps/metamaps/R/dbquery.R")

idk=rbindlist(lapply(unique(summaryTaxa1$originID), function(taxa){
dbquery(
query="START taxon=node:ncbitaxid(taxid={taxa}) RETURN taxon.taxid, taxon.name", 
params= list(taxa=taxa), 
cypherurl = "http://metamaps.scelse.nus.edu.sg:7474/db/data/cypher")
}))

summaryTaxa2=merge(summaryTaxa1, idk, by.x="originID","X1")
merge(summaryTaxa2, top500, by.x="X2",by.y=")

p0=	ggplot(summaryTaxa2, aes(x=reorder(X2,sizePerChild,median), y=log10(sizePerChild)))+
	geom_boxplot()+
	geom_point(alpha=0.3)+
	theme(axis.text.x = element_text(angle=90, size=5, hjust = 1))+
	xlab("Genus")+ylab("# refseq sequences (log10)")
ggsave(plot=p0, file="out/sim.0109.out.pdf",w=20)

#95 taxa recovered
set.seed(2)
new.data=do.call(rbind,
lapply(unique(combined$originID), function(x) { 
    data2 = subset(combined, originID== x)
    data2[sample(nrow(data2),1),]
    }))

write.table(new.data,file="out/sim.0109.out.txt",row.names=F,quote=F,sep="\t")





#-- only 95 have a refseq sequence
#parentchild2=unique(parentchild[,c("originID","originName")])
#parentchild2[parentchild2$originName %in% parentchild2$originName[duplicated(parentchild2$originName)],]
#listofweirdos=unique(parentchild$originID[!parentchild$originID %in% combined$originID])
#[1]  5884  5890  5782 79255 55087 33057 12234

#length(unique(parentchild$originID))
#102

#Reason 
##	1.eukaryota
##	2.duplicate names
#c(1386,55087,79255,2053)

#idk2=rbindlist(lapply(listofweirdos,function(x) { 
#dbquery(
###query = "START taxon=node:ncbitaxid(taxid={taxa}) RETURN taxon.taxid,taxon.name,head(labels(taxon))",
##	params = list(taxa= x), 
#    	query='start basetaxa=node:ncbitaxid(taxid={taxids}) match path=basetaxa-[:childof*]->(superkingdom:superkingdom) return extract(n in nodes(path)| n.name) as name,
# extract(n in nodes(path)| head(labels(n))) as rank',
# 	params = list(taxids=x),	
#	cypherurl = "http://metamaps.scelse.nus.edu.sg:7474/db/data/cypher"
#	)
#}))
