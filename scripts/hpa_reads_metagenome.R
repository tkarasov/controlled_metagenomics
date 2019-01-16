#This script takes the output for HPA and graphs the read progression
#samples that need to be corrected in name before processing:
#D.5.I.5_RunId0091_LaneId1
#control_RunId0091_LaneId1 removed

library(reshape2)
library('dplyr')
library('tidyr')
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
'%!in%' <- function(x,y)!('%in%'(x,y))
mycolors=c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "firebrick", "khaki2", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue","royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")

meta=read.table("~/work_main/abt6_projects9/metagenomic_controlled/data/processed_reads/hpa_infections/meta_family_corrected_per_plant.csv", sep=",", header=T, row.names = 1)
meta_hpa=melt(meta[c("Peronosporaceae"),], value.name="load", variable.name="Genotype")
#meta_pseud$Genotype=rownames(meta_pseud)

#exclude control
meta_hpa=meta_hpa[-c(meta_hpa$Genotype=="control_RunId0091_LaneId1"),]

meta_hpa$Day=as.numeric(as.character(gsub("D", "",sapply(strsplit(as.character(meta_hpa$Genotype), ".", fixed=TRUE), `[`, 1))))
meta_hpa$Replicate=sapply(strsplit(as.character(meta_hpa$Genotype), ".", fixed=TRUE), `[`,3)
meta_hpa$Genotype=sapply(strsplit(as.character(meta_hpa$Genotype), ".", fixed=TRUE), `[`,2)
#meta_pseud$Day=sapply(strsplit(sapply(strsplit(meta_pseud$day_ori, ""), `[`, 1), "_R"), `[`,2)
load_plot=meta_hpa%>% group_by(Day, Genotype) %>% summarise(mean = mean(load, na.rm=T), se = sd(load, na.rm=T) / length(load))

plot_hpa_time <- ggplot(load_plot, aes(x=Day, y=mean, group=Genotype, colour=Genotype))+geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  geom_point()+theme_bw()+geom_line()+theme_bw()+ scale_color_manual(values = wes_palette(n=2, name="GrandBudapest1"))+ylab("Hpa coverage/Plant coverage")

pdf("~/Dropbox/controlled_metagenomics/results_figures/hpa_over_time.pdf",  useDingbats=FALSE, width=6, height=4)
plot_hpa_time + scale_x_continuous(breaks = c(0:11))
dev.off()
#Now for plotting the stacked barplots
#Choose ten most abundant taxa
meta=within(meta, rm("control_RunId0091_LaneId1"))
top10=names(sort(rowSums(meta, na.rm=TRUE), decreasing=TRUE )[1:10])
meta_microbiome=meta[top10,]
rest=colSums(meta[rownames(meta) %!in% top10, ], na.rm=TRUE)
meta_microbiome=rbind(meta_microbiome, rest)
rownames(meta_microbiome)[11]="Rest"
meta_microbiome$Family=rownames(meta_microbiome)
microb_melt=melt(meta_microbiome)
microb_melt$Day=(as.character(gsub("D", "",sapply(strsplit(as.character(microb_melt$variable), ".", fixed=TRUE), `[`, 1))))
microb_melt$Replicate=sapply(strsplit(as.character(microb_melt$variable), ".", fixed=TRUE), `[`,3)
microb_melt$Replicate=sapply(strsplit(as.character(microb_melt$Replicate), "_", fixed=TRUE), `[`,1)
microb_melt$Genotype=sapply(strsplit(as.character(microb_melt$variable), ".", fixed=TRUE), `[`,2)
microb_melt$Day=replace(microb_melt$Day, microb_melt$Day=="11", "664")
microb_melt=microb_melt[order(microb_melt$Day),]
microb_melt=microb_melt[order(microb_melt$Genotype),]
#microb_melt$Day=replace(microb_melt$Day, microb_melt$Day=="11", as.factor("264"))
microb_melt$Gen_Day=paste(microb_melt$Genotype, microb_melt$Day, microb_melt$Replicate, sep=".")


#microb_melt$Day2=relevel(microb_melt$Day, "0")
#microb_melt$Genotype2=relevel(as.factor(microb_melt$Genotype), "C")
#p <- ggplot(data=microb_melt, aes(x=Gen_Day, y=value, fill=Family))
#p + geom_bar(aes(), stat="identity", position="stack") +
#  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +theme(legend.position="bottom", panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(fill=guide_legend(nrow=5))

pdf("~/Dropbox/controlled_metagenomics/results_figures/hpa_time_control.pdf")
microb_melt_C=microb_melt[microb_melt$Genotype=="C",]
p2<- ggplot(data=microb_melt_C, aes(x=Gen_Day, y=value, fill=Family))
p2 + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values =mycolors) +
  theme(legend.position="bottom", panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(fill=guide_legend(nrow=5))+ylim(c(0,0.5))+ylab("Coverage Peronospora/Coverage Plant")
dev.off()
 
pdf("~/Dropbox/controlled_metagenomics/results_figures/hpa_time_treated.pdf")
microb_melt_I=microb_melt[microb_melt$Genotype=="I",]
p1<- ggplot(data=microb_melt_I, aes(x=Gen_Day, y=value, fill=Family))
p1 + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = mycolors) +
  theme(legend.position="bottom", panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(fill=guide_legend(nrow=5))+ylim(c(0,0.5)) +xlab("Day")+ylab("Coverage Peronospora/Coverage Plant")
dev.off()

#make a column for Pseudomonas
per=microb_melt[microb_melt$Family=="Peronosporaceae",]
per$per_val=per$value
per=subset(per,select=-c(Family, value))
meb=merge(microb_melt, per)
meb=meb[which(meb$Genotype!="control"),]
#meb=meb[which(meb$Genotype!="C"),]
lm_Pseud<-lm(log10(value+0.0001)~log10(per_val+0.0001), data=meb[meb$Family=="Pseudomonadaceae",])
lm_Per<-lm(log10(value+0.0001)~log10(per_val+0.0001), data=meb[meb$Family=="Peronosporaceae",])
lm_Ent<-lm(log10(value+0.0001)~log10(per_val+0.0001), data=meb[meb$Family=="Enterobacteriaceae",])
lm_Sphing<-lm(log10(value+0.0001)~log10(per_val+0.0001), data=meb[meb$Family=="Sphingomonadaceae",])
lm_Com<-lm(log10(value+0.0001)~log10(per_val+0.0001), data=meb[meb$Family=="Comamonadaceae",])
lm_Xan<-lm(log10(value+0.0001)~log10(per_val+0.0001), data=meb[meb$Family=="Xanthomonadaceae",])
lm_Rhiz<-lm(log10(value+0.0001)~log10(per_val+0.0001), data=meb[meb$Family=="Rhizobiaceae",])
lm_Morax<-lm(log10(value+0.0001)~log10(per_val+0.0001), data=meb[meb$Family=="Moraxellaceae",])
lm_Brady<-lm(log10(value+0.0001)~log10(per_val+0.0001), data=meb[meb$Family=="Bradyrhizobiaceae",])
lm_Burk<-lm(log10(value+0.0001)~log10(per_val+0.0001), data=meb[meb$Family=="Burkholderiaceae",])
lm_Rest<-lm(log10(value+0.0001)~log10(per_val+0.0001), data=meb[meb$Family=="Rest",])


df=data.frame()
p=ggplot(df) +xlim(-3,1)+ylim(-3,1)

family_vec=list( lm_Brady, lm_Burk, lm_Com, lm_Ent, lm_Morax, lm_Per, lm_Pseud, lm_Rest, lm_Rhiz, lm_Sphing, lm_Xan)

for(i in 1:11){
  family=family_vec[i][[1]]
  slope=as.numeric(as.character(family$coefficients[2]))
  intercept=as.numeric(as.character(family$coefficients[1]))
  if(summary(family)$coefficients[,4][2]<0.005 && i!=6){
    p=p+geom_abline(slope=slope, intercept=intercept, color=mycolors[i]) 
    print(i)
  }
}
pdf("~/Dropbox/controlled_metagenomics/results_figures/effect_hpa_other.pdf")
p+panel_border(colour = "Black",size=1)

dev.off()


#PCA on leftover microbiome

# red is 
red=t((subset(meta_microbiome, select=-c(Family))))#[-c(1),])
#get mean
red_mean=mean(red)


prin_comp=prcomp(red, scale.=T, center=T)
biplot(prin_comp)
ggbiplot(prin_comp)
