#!/usr/bin/env Rscript

library(ggplot2)

data = read.table("out/sim.0302.out.txt",h=T,sep="\t")
completeGenomes = read.table("out/sim.0101.out.txt",sep="\t",h=T)$taxid

data$taxid = do.call(c,lapply(data$ID, function(x) { 
   unlist(strsplit(x=as.character(x), split="\\|"))[2] 
}    ))
data$Genome = 0 
data$Genome[which(data$taxid %in% completeGenomes)] = 1

p0 = ggplot(data, aes(x=as.factor('Chosen'), y=log10(LENGTH)))+
geom_boxplot()+
geom_jitter(aes(color=as.factor(Genome)))+
ylab("Length (log10)")+xlab("Chosen genera")+
ggtitle("Distribution of chosen genome/sequence sizes")+
scale_color_discrete(name="Complete Genome Status")
ggsave(plot=p0, file="out/sim.0303.out.pdf")
save(data, file="out/sim.0303.out.rda")
