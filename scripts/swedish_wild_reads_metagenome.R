library(reshape2)
library('dplyr')
library('tidyr')
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
require(gridExtra)
library(intrval)

meta=read.table("~/work_main/abt6_projects9/metagenomic_controlled/data/processed_reads/swedish_samples/meta_family_corrected_per_plant.csv", sep=",", header=T, row.names = 1)

top10_swed=names(sort(rowSums(meta, na.rm=TRUE), decreasing=TRUE )[1:10])
meta_microbiome_swed=meta[top10_swed,]
rest_swed=colSums(meta[rownames(meta) %ni% top10, ], na.rm=TRUE)
meta_microbiome_swed=rbind(meta_microbiome_swed, rest_swed)
rownames(meta_microbiome_swed)[11]="Rest"
tot_swed=colSums(meta_microbiome_swed)
meta_microbiome_swed$Family=rownames(meta_microbiome_swed)
microb_melt_swed=melt(meta_microbiome_swed, id=c("Family"))
microb_melt_swed$Load=as.factor(-1*tot_swed[microb_melt_swed$variable])

#microb_melt$Day2=relevel(microb_melt$Day, "0")
#microb_melt$Genotype2=relevel(as.factor(microb_melt$Genotype), "C")

manual_colors=c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
                "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")


pdf("~/Dropbox/controlled_metagenomics/results_figures/swedish_metagenome.pdf")
p <- ggplot(data=microb_melt_swed, aes(x=Load, y=value, fill=Family))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = manual_colors) +
  theme(legend.position="bottom", panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank()) + guides(fill=guide_legend(nrow=5)) +
  xlab("Plant Individuals") + ylab("Microbial Cov./Plant Cov.")+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

dev.off()


#Day3 avrB log10     -0.276 
#Day3 C log10       -1.65  
#Day3 EV  log10      1.54  

#Hpa
#Day11 0.0714 or -1.146302 (log10)
#Day5 0.00303 or -2.518557 (log10)
#How much Pseudomonas?
pdf("~/Dropbox/controlled_metagenomics/results_figures/sweden_pseud_hpa.pdf")
families=data.frame(t(meta))
pseud=ggplot(data=families, aes(log10(Pseudomonadaceae))) + geom_histogram(fill="GRAY", colour = "BLACK")+theme_bw()+xlab(expression(log[10]~("Pseudomonadaceae coverage"/"Plant coverage")))+geom_vline(aes(xintercept=-0.276, color="Resistant"))+geom_vline(aes(xintercept=1.54, color="Susceptible"))
hpa=ggplot(data=families, aes(log10(Peronosporaceae))) + geom_histogram(fill="GRAY", colour = "BLACK")+theme_bw()+xlab(expression(log[10]~("Peronosporaceae coverage"/"Plant coverage")))+geom_vline(aes(xintercept=-1.146302, color="Day 11"))+geom_vline(aes(xintercept=-2.518557 , color="Day 5"))
grid.arrange(pseud, hpa, ncol=1)
dev.off()
