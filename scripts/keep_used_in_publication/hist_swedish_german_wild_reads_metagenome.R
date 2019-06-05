library(reshape2)
library('dplyr')
library('tidyr')
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(matrixStats)
require(gridExtra)
library(cowplot)
#library(moments)
library(intrval)
library("Hmisc")

manual_colors=c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")


#The goal of this script is to generate figures side-by-side that will take the metagenome corrected tables for germany and sweden and graph the bar plots next to one another

meta_table<-function(level){
  concat_name=paste(level,".csv",sep="")
  meta_sweden=read.table(paste("~/work_main/abt6_projects9/metagenomic_controlled/data/processed_reads/swedish_samples/meta_family_corrected_per_plant_", concat_name, sep=""), sep=",",header=T, row.names = 1)
  meta_german1=read.table(paste("~/work_main/abt6_projects9/metagenomic_controlled/data/processed_reads/german_samples/meta_family_corrected_per_plant_", concat_name, sep=""), sep=",",header=T, row.names = 1)
  meta_german=meta_german1[,keep$metagenome_identifier]
  top_swed=sort(rowSums(meta_sweden, na.rm=TRUE), decreasing=TRUE )[1:5]
  top_germ=sort(rowSums(meta_german, na.rm=TRUE), decreasing=TRUE )[1:5]
  together=unique(c(names(top_germ), names(top_swed)))
  germany_processed=process_microb_melt(meta_german, together)
  germany_processed$value=as.numeric(as.character(germany_processed$value))
  sweden_processed=process_microb_melt(meta_sweden,together)
  sweden_processed$value=as.numeric(as.character(sweden_processed$value))
  out<-list()
  out$germany=germany_processed
  out$sweden=sweden_processed
  return(out)
}

process_microb_melt<-function(meta, together){
  meta_microbiome=meta[together,]
  rownames(meta_microbiome)=c(together)
  rest=colSums(meta[rownames(meta) %ni% together, ], na.rm=TRUE)
  meta_microbiome=rbind(meta_microbiome, rest)
  rownames(meta_microbiome)[dim(meta_microbiome)[1]]="Rest"
  tot=colSums(meta_microbiome, na.rm=TRUE)
  meta_microbiome$Family=rownames(meta_microbiome)
  microb_melt=melt(meta_microbiome, id=c("Family"))
  #microb_melt[is.na(microb_melt$value),]$value=0
  microb_melt$Load=as.factor(-1*tot[microb_melt$variable])
  return(microb_melt)
}

plot_side_by_side<-function(out){
  max_y_axis=max(max(out$germany[,'value']), max(out$sweden[,'value']))
  
  b_g=ggplot(data=out$germany, aes(x=Load, y=value, fill=Family))+ geom_bar(aes(), stat="identity", position="stack") +
    scale_fill_manual(values = manual_colors) +
    theme(legend.position="bottom", panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank()) + 
    guides(fill=guide_legend(nrow=5)) +
    xlab("Plant Individuals") + ylab("Microbial Cov./Plant Cov.") +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_y_axis)) +
    annotate("text",  x=Inf, y = Inf, label = "Germany", vjust=2, hjust=2)
    
    s_g=ggplot(data=out$sweden, aes(x=Load, y=value, fill=Family))+ geom_bar(aes(), stat="identity", position="stack") +
      scale_fill_manual(values = manual_colors) +
      theme(legend.position="bottom", panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank()) + 
      guides(fill=guide_legend(nrow=5)) +
      xlab("Plant Individuals") + ylab("Microbial Cov./Plant Cov.") +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
      scale_y_continuous(expand = c(0, 0), limits = c(0, max_y_axis)) +
      annotate("text",  x=Inf, y = Inf, label = "Sweden", vjust=2, hjust=2) 
    
    grobs <- ggplotGrob(s_g)$grobs
    legend<-grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
    pgrid<-plot_grid(b_g+theme(legend.position="none"), s_g+theme(legend.position="none")+
                       theme(axis.title.y=element_blank()))
    p <- plot_grid(pgrid, legend, nrow = 2, ncol=1, rel_heights = c(1, .3))
    
    p
}

#first thing is to remove the offending columns from the german dataset
german_samples=read.table("~/Dropbox/germany_pathogen_collections/sample_data/plate_sample_locations/sample_infoFinal_2018.txt", sep="\t", header=T)
keep=german_samples[-which(duplicated(german_samples$uniqueID)),]


#generate_tables
out_bacteria=meta_table("bacteria")
p_bac=plot_side_by_side(out_bacteria)

out_fungi=meta_table("fungi")
p_fungi=plot_side_by_side(out_fungi)

out_oom=meta_table("oomycete")
p_oom=plot_side_by_side(out_oom)

p_tot=plot_grid(p_bac, p_oom, p_fungi, nrow)

pdf("~/Dropbox/controlled_metagenomics/results_figures/oomycete_germany_sweden.pdf")
p_oom
dev.off()

pdf("~/Dropbox/controlled_metagenomics/results_figures/bacteria_germany_sweden.pdf")
p_bac
dev.off()

pdf("~/Dropbox/controlled_metagenomics/results_figures/fungi_germany_sweden.pdf")
p_fungi
dev.off()


hm=read.csv("/Users/tkarasov/work_main/abt6_projects9/metagenomic_controlled/data/processed_reads/german_samples/meta_family_corrected_per_plant.csv", sep=",",header=T, row.names = 1)
hm_pseud=as.numeric(as.character((hm["Pseudomonadaceae",])))
germ_pseud=as.numeric(as.character(meta_german["Pseudomonadaceae",]))

pdf("~/Dropbox/controlled_metagenomics/results_figures/german_metagenome.pdf")
p <- ggplot(data=microb_melt, aes(x=Load, y=value, fill=Family))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = manual_colors) +
  theme(legend.position="bottom", panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank()) + 
  guides(fill=guide_legend(nrow=5)) +
  xlab("Plant Individuals") + ylab("Microbial Cov./Plant Cov.") +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 225)) 

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
