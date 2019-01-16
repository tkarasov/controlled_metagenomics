library(reshape2)
library('dplyr')
library('tidyr')
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(matixStats)
require(gridExtra)
library(moments)

meta=read.table("~/work_main/abt6_projects9/metagenomic_controlled/data/processed_reads/german_samples/meta_family_corrected_per_plant.csv", sep=",", header=T, row.names = 1)

top10=names(sort(rowSums(meta, na.rm=TRUE), decreasing=TRUE )[1:10])
meta_microbiome=meta[top10,]
rest=colSums(meta[rownames(meta) %!in% top10, ], na.rm=TRUE)
meta_microbiome=rbind(meta_microbiome, rest)
rownames(meta_microbiome)[11]="Rest"
tot=colSums(meta_microbiome)
meta_microbiome$Family=rownames(meta_microbiome)
microb_melt=melt(meta_microbiome, id=c("Family"))
microb_melt$Load=as.factor(-1*tot[microb_melt$variable])

#microb_melt$Day2=relevel(microb_melt$Day, "0")
#microb_melt$Genotype2=relevel(as.factor(microb_melt$Genotype), "C")
pdf("~/Dropbox/controlled_metagenomics/results_figures/german_metagenome.pdf")
p <- ggplot(data=microb_melt, aes(x=Load, y=value, fill=Family))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
                               "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  theme(legend.position="bottom", panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank()) + guides(fill=guide_legend(nrow=5)) +
  xlab("Plant Individuals") + ylab("Microbial Cov./Plant Cov.")

dev.off()

pdf("~/Dropbox/controlled_metagenomics/results_figures/germany_pseud_hpa.pdf")
families=data.frame(t(meta))
pseud=ggplot(data=families, aes(log10(Pseudomonadaceae))) + geom_histogram(fill="GRAY", colour = "BLACK")+theme_bw()+xlab(expression(log[10]~("Pseudomonadaceae coverage"/"Plant coverage")))+geom_vline(aes(xintercept=-0.276, color="Resistant"))+geom_vline(aes(xintercept=1.54, color="Susceptible"))
hpa=ggplot(data=families, aes(log10(Peronosporaceae))) + geom_histogram(fill="GRAY", colour = "BLACK")+theme_bw()+xlab(expression(log[10]~("Peronosporaceae coverage"/"Plant coverage")))+geom_vline(aes(xintercept=-1.146302, color="Day 11"))+geom_vline(aes(xintercept=-2.518557 , color="Day 5"))
grid.arrange(pseud, hpa, ncol=1)
dev.off()

#Now plot variance and mean statistics
MatVar <- function(x, dim = 1, ...) {
  if(dim == 1){
    rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
  } else if (dim == 2) {
    rowSums((t(x) - colMeans(x, ...))^2, ...)/(dim(x)[1] - 1)
  } else stop("Please enter valid dimension")
}

sd_meta=sqrt(MatVar(meta))
skew=skewness(t(meta), na.rm=T)
kurt=kurtosis(t(meta),na.rm=T)
mean_meta=rowSums(meta)/dim(meta)[1]
together=cbind(mean_meta, sd_meta)
max_t=apply(meta, 1, max)
together=cbind(together, max_t)
together=cbind(together, skew)
together=cbind(together, kurt)
together=as.data.frame(together[order(-mean_meta),])
var_plot<-ggplot(data=together, aes(x=c(1:dim(together)[1]),y=mean_meta))

var_plot + geom_point() + geom_errorbar(ymin=mean_meta-sd_meta, ymax=mean_meta+sd_meta, width=.1)

quant=rowQuantiles(as.matrix(meta))
