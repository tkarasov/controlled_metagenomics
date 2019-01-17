library(reshape2)
library('dplyr')
library('tidyr')
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
require(gridExtra)
library(cowplot)
library(lmtest)
'%!in%' <- function(x,y)!('%in%'(x,y))
mycolors=c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "firebrick", "khaki2", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue","royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")

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

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#meta=read.table("~/work_main/abt6_projects9/metagenomic_controlled/data/processed_reads/dc3000_infections/meta_family_corrected_per_plant.csv", sep=",", header=T, row.names = 1)
meta=read.table("/ebio/abt6_projects9/metagenomic_controlled/data/processed_reads/dc3000_infections/meta_family_corrected_per_plant.csv", sep=",", header=T, row.names = 1)


top10=names(sort(rowSums(meta, na.rm=TRUE), decreasing=TRUE )[1:10])
meta_microbiome=meta[top10,]
rest=colSums(meta[rownames(meta) %!in% top10, ], na.rm=TRUE)
meta_microbiome=rbind(meta_microbiome, rest)
rownames(meta_microbiome)[11]="Rest"
tot=colSums(meta_microbiome)
meta_microbiome$Family=rownames(meta_microbiome)
microb_melt=melt(meta_microbiome, id=c("Family"))
microb_melt$day_ori=sapply(strsplit(as.character(microb_melt$variable), "_"), `[`, 2)
microb_melt$Genotype=sapply(strsplit(as.character(microb_melt$variable), "_"), `[`,1)
microb_melt$Day=sapply(strsplit(microb_melt$day_ori, "\\."), `[`, 1)
microb_melt$Day=gsub("R","", microb_melt$Day)
microb_melt$Day=sapply(microb_melt$Day, replace_rr)
#microb_melt$Genotype=gsub("avrB", "Nope", microb_melt$Genotype)
microb_melt$Replicate=sapply(strsplit(microb_melt$day_ori, "\\."), `[`, 2)
microb_melt$Replicate=gsub( "_","", microb_melt$Replicate)
microb_melt$combined=paste(microb_melt$Day, microb_melt$Replicate, sep="_")

#ONLY EV
EV=microb_melt[microb_melt$Genotype=="EV",]
EV=EV[which(is.na(EV$Replicate)==FALSE),]
EV_plot <- ggplot(data=EV, aes(x=combined, y=log10(value+1), fill=Family))

########
avrB=microb_melt[microb_melt$Genotype=="avrB",]
avrB=avrB[which(is.na(avrB$Replicate)==FALSE),]
avrB_plot <- ggplot(data=avrB, aes(x=combined, y=log10(value+1), fill=Family))

########
control=microb_melt[microb_melt$Genotype=="C",]
control=control[which(is.na(control$Replicate)==FALSE),]
control_plot <- ggplot(data=control, aes(x=combined, y=log10(value+1), fill=Family))




###
avrB_final= avrB_plot + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = mycolors) +
  theme(legend.position="none", panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(fill=guide_legend(nrow=5)) +
  xlab("") +  ylab(expression(log[10]~("Microbe cov."/"Plant cov.")))+ ylim(c(0,3))

control_final= control_plot + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = mycolors) +
  theme(legend.position="none", panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(fill=guide_legend(nrow=5)) +
  xlab("") +  ylim(c(0,3)) + ylab("")

with_legend<- control_plot + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = mycolors) +
  theme(legend.position="bottom", panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(fill=guide_legend(nrow=5)) +
  xlab("Plant Individuals") +  ylab(expression(log[10]~("Microbe cov."/"Plant cov."))) + ylim(c(0,3))

EV_final= EV_plot + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = mycolors) +
  theme(legend.position="none", panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(fill=guide_legend(nrow=5)) +
  xlab("Day") + ylab("")+ ylim(c(0,3))


leg<-g_legend(with_legend)

pdf("~/Dropbox/controlled_metagenomics/results_figures/dc3000_temporal_metagenome_Triplot.pdf", width=8, height=12)
plot_grid(control_final, avrB_final, EV_final, leg, rel_heights=c(1/3.5,1/3.5,1/3.5, .5/3.5),  nrow=4)
dev.off()


if(FALSE) {
pdf("~/Dropbox/controlled_metagenomics/results_figures/dc3000_metagenome.pdf")
p <- ggplot(data=microb_melt, aes(x=combined, y=value, fill=Family))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = mycolors) +
theme(legend.position="bottom", panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(fill=guide_legend(nrow=5)) +
  xlab("Plant Individuals") + ylab("Microbial Cov./Plant Cov.")

dev.off()
}


#now regression modeling
#beyond day 2 the infection seems to plateau
microb_melt$infect_date=microb_melt$Day>2
microb_melt$Genotype=as.factor(microb_melt$Genotype)
microb_melt$Genotype<-relevel(microb_melt$Genotype, ref=2)

#make a column for Pseudomonas
Pseud=microb_melt[microb_melt$Family=="Pseudomonadaceae",]
Pseud$pseud=Pseud$value
Pseud=subset(Pseud, select=-c(Family, value))
meb=merge(microb_melt, Pseud)
meb=meb[which(meb$Genotype!="control"),]
meb=meb[which(meb$Genotype=="EV"),]
#meb=meb[which(meb$Genotype!="C"),]
lm_Pseud<-lm(log10(value)~log10(pseud), data=meb[meb$Family=="Pseudomonadaceae",])
lm_Ent<-lm(log10(value)~log10(pseud), data=meb[meb$Family=="Enterobacteriaceae",])
lm_Caul<-lm(log10(value)~log10(pseud), data=meb[meb$Family=="Caulobacteraceae",])
lm_Sphing<-lm(log10(value)~log10(pseud), data=meb[meb$Family=="Sphingomonadaceae",])
lm_Alcal<-lm(log10(value)~log10(pseud), data=meb[meb$Family=="Alcaligenaceae",])
lm_Xan<-lm(log10(value)~log10(pseud), data=meb[meb$Family=="Xanthomonadaceae",])
lm_Rhiz<-lm(log10(value)~log10(pseud), data=meb[meb$Family=="Rhizobiaceae",])
lm_Morax<-lm(log10(value)~log10(pseud), data=meb[meb$Family=="Moraxellaceae",])
lm_Brady<-lm(log10(value)~log10(pseud), data=meb[meb$Family=="Bradyrhizobiaceae",])
lm_Bruc<-lm(log10(value)~log10(pseud), data=meb[meb$Family=="Brucellaceae",])
lm_Rest<-lm(log10(value)~log10(pseud), data=meb[meb$Family=="Rest",])

#now fit second regression model with squared residuals of first. This is how one would test for heteroskedasticity but I dont see it. Hm
conditonal_variance<-lmtest::bptest(lm_Ent)


###
df=data.frame()
p=ggplot(df) +xlim(-2,3)+ylim(-2,0.0)

family_vec=list(lm_Alcal, lm_Brady, lm_Bruc, lm_Caul, lm_Ent, lm_Morax, lm_Pseud, lm_Rest, lm_Rhiz, lm_Sphing, lm_Xan)

for(i in 1:11){
  family=family_vec[i][[1]]
  slope=as.numeric(as.character(family$coefficients[2]))
  intercept=as.numeric(as.character(family$coefficients[1]))
  if(summary(family)$coefficients[,4][2]<0.005 && i!=7){
    p=p+geom_abline(slope=slope, intercept=intercept, color=mycolors[i]) 
    print(i)
  }
  }
pdf("~/Dropbox/controlled_metagenomics/results_figures/effect_dc3000_other.pdf")
p+panel_border(colour = "Black",size=1)

dev.off()


plot_regress=function(name_family){
  nam = paste("p_",name_family, sep="")
  #assign(nam, ggplot(data=meb[meb$Family==name_family,], aes(x=per_val, y=value)) +
  #         geom_point() + geom_smooth(method=lm, se=TRUE, col="BLACK") + theme_bw())
  my_subset=meb[meb$Family==name_family,]
  regress = lm(log10(value+0.000001)~log10(pseud), data=my_subset)
  is_sig=summary(regress)$coefficients[,4][2]
  sig=FALSE
  basic= ggplot(data=my_subset, aes(x=log10(pseud), y=log10(value))) +
    geom_point(cex=0.2) + geom_smooth(method=lm, se=TRUE, col="BLACK") + theme_bw()
  fit = lm(log10(value+0.000001)~log10(pseud), data=my_subset)
  rsquared=round(signif(summary(fit)$adj.r.squared, 5), digits=3)
  rsq_nam=paste(" R2",rsquared,sep="=")
  name_family_amend = paste(name_family, rsq_nam, sep=",")#paste(name_family, rsquared, sep=rsq_nam)
  
  if(is_sig<0.005){
    p= basic + annotate("text", label=paste(name_family_amend, "*", sep=""), x=0, y=1, cex=2, color="RED") +
      xlim(-2.5,2) + ylim(-5,2) +xlab("") + ylab("")#+ xlab(expression(log[10]~("Hpa Load"))) + ylab(expression(log[10]~("Other Load")))
    return(p)}
  else{
    p = basic + annotate("text", label=name_family_amend, x=0, y=1, cex=2, color="RED") +
      xlim(-2.5,2) + ylim(-5,2) +xlab("") +ylab("") # + xlab(expression(log[10]~("Hpa Load"))) + ylab(expression(log[10]~("Other Load"))) 
    return(p)
  }
}
families=c("Pseudomonadaceae","Enterobacteriaceae","Sphingomonadaceae","Caulobacteraceae","Xanthomonadaceae","Rhizobiaceae","Moraxellaceae","Bradyrhizobiaceae","Alcaligenaceae", "Brucellaceae","Rest")
pdf("/ebio/abt6_projects9/metagenomic_controlled/code/controlled_metagenomics_git/data/dc3000_other_taxa_correlations.pdf")
i=1
plots=rep(0, length(families))
for(fam in families){
  print(fam)
  assign(paste("p_", i, sep=""), plot_regress(fam))
  i=i+1
}
grid.arrange(p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9,p_10, bottom=textGrob(expression(log[10]~"(Pseud. Cov./Plant Cov.)")), left=textGrob(expression(log[10]~"(Focal Cov./Plant Cov.)"), rot=90))
dev.off()

#Rank correlations vs. other for 
pdf("/ebio/abt6_projects9/metagenomic_controlled/code/controlled_metagenomics_git/data/hist_dc3000_hpa.pdf)
EV_only = grep("EV_", colnames(meta))
meta_EV=t(meta[,EV_only])
hpa=read.table("/ebio/abt6_projects9/metagenomic_controlled/data/processed_reads/hpa_infections//meta_family_corrected_per_plant.csv", sep=",", header=T, row.names = 1)
hpa_I_only=grep(".I.", colnames(hpa))
hpa_EV=t(hpa[,hpa_I_only])
corr_pseud=cor(meta_EV, method='spearman')
corr_hpa=cor(hpa_EV, method='spearman')
per_fin=t(melt(corr_hpa[,"Peronosporaceae"], value.name="per"))
pseud=data.frame(t(melt(corr_pseud[,"Pseudomonadaceae"],value.name="pseud")))
per_pseud=data.frame(t(rbind(per_fin[, intersect(colnames(per_fin), colnames(pseud))], pseud[,intersect(colnames(per_fin), colnames(pseud))])))
colnames(per_pseud)=c("Peronosporaceae","Pseudomonadaceae")
per_pseud=per_pseud[which(per_pseud$Peronosporaceae!="NA"),]
per_pseud=per_pseud[which(per_pseud$Pseudomonadaceae!="NA"),]
per_pseud=per_pseud[-c(which(rownames(per_pseud)=="Peronosporaceae")),]
per_pseud=per_pseud[-c(which(rownames(per_pseud)=="Pseudomonadaceae")),]
m_ps=melt(per_pseud)


plot_hist = ggplot(data=m_ps, aes(x=value, fill=variable, stat(density)))

plot_hist+geom_histogram(alpha=1)+scale_fill_brewer(palette="Set1") + xlab("Spearman rank correlation, R") + xlim(-1,1)
dev.off()






red=t((subset(meta_microbiome, select=-c(Family, control_1, control_2, control_2, control_3, control_4))))#[-c(1),])
#get mean
red_mean=mean(red)


prin_comp=prcomp(red, scale.=T, center=T)
biplot(prin_comp)
ggbiplot(prin_comp)
