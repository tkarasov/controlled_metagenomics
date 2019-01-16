library(reshape2)
library('dplyr')
library('tidyr')
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
require(gridExtra)

replace_rr <- function(rr){
  if(is.na(rr)) return(NA)
  else if(rr=="I") return(0)
  else if(rr=="II") return(1)
  else if(rr=="III") return(2)
  else if(rr=="IV") return(3)
  else if(rr=="V") return(4)
  else if(rr=="VI") return(6)
  else return(NA)
}



#Combining different metrics
colony_count=read.table("~/Dropbox/controlled_metagenomics/infection_results/dc3000_colony_count.txt", sep="\t")

qPCR = read.table("~/Dropbox/controlled_metagenomics/qPCR/qPCR_delta.txt", sep="\t")
colnames(qPCR)[1:3]=c("Genotype", "Day", "Replicate")

qPCR=qPCR[, c("Day", "Replicate", "Genotype",  "ratio")]

colony_count$Genotype=gsub(" ", "",colony_count$Genotype)


combined = dplyr::full_join(colony_count, qPCR)
ggplot(data=combined, aes(x=log10(FINAL+.00001), y=log10(ratio+0.00001))) +  geom_point() + xlab(expression(log[10]~(cfu/mL))) + ylab(expression(log[10]~("ng bacteria"/"ng plant"))) + geom_smooth(method=lm, se = FALSE, col="Gray")
       
cor.test(log10(combined$ratio+0.000001), log10(combined$FINAL+.000001))

#metagenomics data

meta=read.table("~/work_main/abt6_projects9/metagenomic_controlled/data/processed_reads/dc3000_infections/meta_family_corrected_per_plant.csv", sep=",", header=T, row.names = 1)
meta_pseud=melt(meta[c("Pseudomonadaceae"),], value.name="load", variable.name="Genotype")
#meta_pseud$Genotype=rownames(meta_pseud)
meta_pseud$day_ori=sapply(strsplit(as.character(meta_pseud$Genotype), "-"), `[`, 1)
meta_pseud$Genotype=sapply(strsplit(as.character(meta_pseud$Genotype), "_"), `[`,1)
meta_pseud$Day=sapply(strsplit(sapply(strsplit(meta_pseud$day_ori, "\\."), `[`, 1), "_R"), `[`,2)
meta_pseud$Replicate=sapply(strsplit(meta_pseud$day_ori, "\\."), `[`, 2)
meta_pseud$Replicate=gsub( "_","", meta_pseud$Replicate)
meta_pseud$Day=sapply(meta_pseud$Day, replace_rr)
all_all=merge(combined,meta_pseud)
all_keep=all_all[which(all_all$ratio!="NA"),]
ggplot(data=all_keep, aes(x=log10(FINAL+1), y=log10(load))) +  geom_point() + xlab(expression(log[10]~(cfu/mL))) + ylab(expression(log[10]~("ng bacteria"/"ng plant"))) + geom_smooth(method=lm, se = FALSE, col="Gray")

load_plot=all_keep%>% group_by(Day, Genotype) %>% summarise(mean = mean(log10(load+0.0001), na.rm=T), se = sd(log10(load+.0001), na.rm=T) / length(load)) 
cfu_plot=all_keep%>% group_by(Day, Genotype) %>% summarise(mean = mean(log10(FINAL+1), na.rm=T), se = sd(log10(FINAL+1), na.rm=T) / length(FINAL)) 

qpcr_plot=all_keep%>% group_by(Day, Genotype) %>% summarise(mean = mean(log10(ratio+0.0001), na.rm=T), se = sd(log10(ratio+0.0001), na.rm=T) / length(ratio)) 


df <- bind_rows(CFU=cfu_plot, qPCR=qpcr_plot, Metagenome=load_plot, .id="group")

pdf("~/Dropbox/controlled_metagenomics/results_figures/tri_comparison.pdf")
ggplot(df, aes(x=Day, y=mean, group=Genotype, colour=Genotype))+geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  geom_point()+theme_bw()+geom_line()+theme_bw()+ 
  facet_wrap(~group , scales = "free_y", ncol=1, strip.position = "right" )             + 
  ylab(NULL) + theme(legend.position = "none") + scale_color_manual(values = wes_palette(n=3, name="GrandBudapest1"))
dev.off()

#plot each against each!
pdf("~/Dropbox/controlled_metagenomics/results_figures/cfu_vs_other.pdf", height=12, width=6)
plot1 <- ggplot(data=all_keep, aes(x=log10(FINAL+1), y=log10(ratio+0.01)))+geom_point()  + xlab(expression(log[10]~(cfu/mL))) + ylab(expression(log[10]~("ng bact."/"ng plant"))) + geom_smooth(method=lm, se = FALSE, col="Gray") + theme_bw()
plot2 <- ggplot(data=all_keep, aes(x=log10(FINAL+1), y=log10(load+0.01)))+geom_point() + xlab(expression(log[10]~(cfu/mL))) + ylab(expression(log[10]~("DC3000 coverage"/"Plant coverage"))) + geom_smooth(method=lm, se = FALSE, col="Gray") +theme_bw()
plot3 <- ggplot(data=all_keep, aes(x=log10(ratio+0.01), y=log10(load+0.01)))+geom_point() + ylab(expression(log[10]~("DC3000 coverage"/"Plant coverage"))) + xlab(expression(log[10]~("ng bact."/"ng plant"))) + geom_smooth(method=lm, se = FALSE, col="Gray") +theme_bw()
grid.arrange(plot1, plot2,plot3, ncol=1)
dev.off()

cor.test(log10(all_keep$ratio+0.01), log10(all_keep$FINAL+1))
#0.612
cor.test(log10(all_keep$load+0.01), log10(all_keep$FINAL+1))
#.767


#now for plotting the metagenome
#microb_melt$Day2=relevel(microb_melt$Day, "0")
#microb_melt$Genotype2=relevel(as.factor(microb_melt$Genotype), "C")
p <- ggplot(data=microb_melt, aes(x=Gen_Day, y=value, fill=Family))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
                               "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  theme(legend.position="bottom", panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(fill=guide_legend(nrow=5))


microb_melt_C=microb_melt[microb_melt$Genotype=="C",]
p2<- ggplot(data=microb_melt_C, aes(x=Gen_Day, y=value, fill=Family))
p2 + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
                               "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  theme(legend.position="bottom", panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(fill=guide_legend(nrow=5))+ylim(c(0,0.5))+ylab("Coverage Peronospora/Coverage Plant")


microb_melt_I=microb_melt[microb_melt$Genotype=="I",]
p1<- ggplot(data=microb_melt_I, aes(x=Gen_Day, y=value, fill=Family))
p1 + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
                               "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  theme(legend.position="bottom", panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(fill=guide_legend(nrow=5))+ylim(c(0,0.5)) +xlab("Day")+ylab("Coverage Peronospora/Coverage Plant")



