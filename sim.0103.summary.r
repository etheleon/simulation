#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

data=read.table("out/sim.0102.out",sep="\t", h=T)
bplength = setNames(read.table("out/sim.0102.out2",h=F, sep="\t"), c("parentGenus", "bplength"))


d2=data %.% 
group_by(parentGenus) %.%
summarise(numSeq=n(),numtaxa=length(unique(taxid)))

d2$ratio=d2$numSeq/d2$numtaxa
d3=merge(d2, bplength, by="parentGenus")


p0 = ggplot(d3, aes(y=numSeq, x=numtaxa))+
geom_point(aes(size=sqrt(bplength)))+
geom_text(data=subset(d3, numSeq > 20000| numtaxa >=50 ), aes(label=parentGenus), color='red')+
ggtitle("Number of Sequences against Number of taxa per Genera")+
xlab("# taxa")+ylab("# sequences")
ggsave(plot=p0, filename="out/sim.0103.summary1.pdf")

pdf("out/sim.0103.summary.pdf",w=10,h=10)
lapply(unique(data$parentGenus), function(x) { 
  df=subset(data, parentGenus == x) 
  df2=df %.% 
  group_by(taxid) %.%
  summarise(numSeq=n())
  qplot(reorder(as.factor(taxid),numSeq), log10(numSeq), data=df2, geom="bar", stat="identity")+
  theme(axis.text.x=element_text(angle=45, hjust=1,vjust=1))+ggtitle(x)
})
dev.off()


##########
#luckydraw
##########

set.seed(5)
chosen = do.call(rbind,lapply(unique(data$parentGenus), function(x) { 
    chosen=sample(size=1, as.character(unique(subset(data, parentGenus == x)$taxid)))
    subset(data, taxid == chosen)
}))
write.table(chosen, file="out/sim.0103.chosen.txt", sep="\t", quote=F, row.names=F)
