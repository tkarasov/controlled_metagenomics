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
library(ape)
library(cowplot)
library(vegan)
library(extrafont)
library(lme4)
#font_import()
loadfonts()

process_gs<-function(g_s, subset=FALSE, classifier=cf, classifier_identity=cfi){
  #function that takes gs
  if(subset==TRUE){
    g_temp<-g_s[which(g_s[[cf]]==cfi),]
  }
  else{
    g_temp<-g_s
  }
  g_s_otus <- g_temp[,which(colnames(g_temp)%ni%c("population", "country", "load"))]
  g_s_load <-g_temp$load[which(apply((g_s_otus), 1, var)!=0)]
  g_s_pop <-g_temp$population[which(apply((g_s_otus), 1, var)!=0)]
  g_s_country <-g_temp$country[which(apply((g_s_otus), 1, var)!=0)]
  #g_s_otus <-g_s[,which(apply((g_s), 1, var)!=0)]
  #first samples that are missing
  g_s_otus <-g_s_otus[which(apply((g_s_otus), 1, var)!=0),]
  #Next taxa that are absent
  g_s_otus <-g_s_otus[which(apply(t(g_s_otus), 2, var)!=0),]
  #remove families that are 
  keep_0.001 <- which(colSums(g_s_otus)/dim(g_s_otus)[1]>0.001)
  g_s_otus_filter <- g_s_otus[,keep_0.001]
  #now make compositional
  g_s_otus_filter <- g_s_otus_filter/rowSums(g_s_otus)
  together<-list(g_s_otus_filter, g_s_country, g_s_pop, g_s_load)
  return(together)
}

plot_pcoa<-function(pcoa_thing){
  ggplot(data=data.frame(pcoa_thing$vectors), aes(x=Axis.1, y=Axis.2)) + 
    geom_point(col= brewer.pal(n = 8, name = "Set1")[as.factor(g_s_country)]) +
    theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab(paste(paste("Axis 1:", (100*round(perc_explained[1],3)), sep=" "), "%", sep="")) +
    ylab(paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))
}



#Now let's look at PCAs and regressions of their comparison

for(concat_name in c("_bacteria.csv", "_fungi.csv", "_oomycete.csv", "_virus.csv")){
  g_s=read.csv(paste("~/Dropbox/controlled_metagenomics/results_figures/sweden_germany_combined", concat_name, sep=''), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  g_s_otus_filter = process_gs(g_s, subset=FALSE, classifier = "NA", classifier_identity = "NA")

g_s_otus <- g_s[,which(colnames(g_s)%ni%c("population", "country", "load"))]
g_s_load <-g_s$load[which(apply((g_s_otus), 1, var)!=0)]
g_s_pop <-g_s$population[which(apply((g_s_otus), 1, var)!=0)]
g_s_country <-g_s$country[which(apply((g_s_otus), 1, var)!=0)]


#first samples that are missing
g_s_otus <-g_s_otus[which(apply((g_s_otus), 1, var)!=0),]

#Next taxa that are absent
g_s_otus <-g_s_otus[which(apply(t(g_s_otus), 2, var)!=0),]

#remove families that are 
keep_0.001 <- which(colSums(g_s_otus)/dim(g_s_otus)[1]>0.001)
g_s_otus_filter <- g_s_otus[,keep_0.001]

#now make compositional
g_s_otus_filter <- g_s_otus_filter/rowSums(g_s_otus)




#####PCA: How do the populations differ in their compositions?
pc<-prcomp((g_s_otus_filter), center = TRUE, scale. = TRUE, retx = TRUE)
ev <- pc$sdev^2
#Jari Oksanen explanation https://stat.ethz.ch/pipermail/r-help/2005-August/076610.html
prop_ev <- ev/sum(ev)
colfunc <- colorRampPalette(c("yellow", "purple"))

#x	if retx is true the value of the rotated data (the centred (and scaled if requested) data multiplied by the rotation matrix) is returned. Hence, cov(x) is the diagonal matrix diag(sdev^2). For the formula method, napredict() is applied to handle the treatment of values omitted by the na.action.

pc1 <- pc$x[,1]
pc2 <- pc$x[,2]
pc3 <- pc$x[,3]
summary(lmer(pc1~(g_s_load)+(1|g_s_pop)))
summary(lm(pc2~g_s_load+g_s_country))
summary(lm(pc3~g_s_load+g_s_country))

pca_vectors = data.frame(pc$x)
ggplot(data=pca_vectors, aes(x=PC1, y=PC2)) + 
  geom_point(col= brewer.pal(n = 8, name = "Set1")[as.factor(g_s_country)]) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(paste(paste("Axis 1:", (100*round(prop_ev[1],3)), sep=" "), "%", sep="")) +
  ylab(paste(paste("Axis 2:", (100*round(prop_ev[2],3)), sep=" "), "%", sep=""))

ggplot(data=pca_vectors, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color=as.numeric(as.character(log10(g_s_load))))) +
  scale_color_gradient2(name=expression(paste('log'[10],"(Tot. Load)", sep="")),  low = "blue", mid = "white",
                        high = "red") +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") + 
  theme(legend.key.width=unit(1,"cm"), axis.title.y = element_blank())
#ylab(paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))


#####BRAY CURTIS: Do these populations differ in who is there?
#now let's think about Bray-Curtis and PCoA
bray_curtis <-as.matrix(vegdist(g_s_otus_filter, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE))

#PCoA on bray-curtis
bc_pcoa=pcoa(bray_curtis, correction = "lingoes")
bip=biplot(bc_pcoa, pch=20)
perc_explained=bc_pcoa$values[,3]

plot_country <- ggplot(data=data.frame(bc_pcoa$vectors), aes(x=Axis.1, y=Axis.2)) + 
  geom_point(col= brewer.pal(n = 8, name = "Set1")[as.factor(g_s_country)]) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(paste(paste("Axis 1:", (100*round(perc_explained[1],3)), sep=" "), "%", sep="")) +
  ylab(paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))

#pdf("~/Dropbox/controlled_metagenomics/results_figures/pcoa_bray_curtis_g_s.pdf")
#plot_pcoa(bc_pcoa)
#dev.off()

####Do PC1 and PC2 distinguish between samples from germany vs sweden
pc1<-bc_pcoa$vectors[,1]
german_pc1<-pc1[g_s_country=="germany"]
sweden_pc1<-pc1[g_s_country=="Sweden"]
wilcox.test(german_pc1, sweden_pc1)
#p=0.0033 W=5367

pc2<-bc_pcoa$vectors[,2]
german_pc2<-pc2[g_s_country=="germany"]
sweden_pc2<-pc2[g_s_country=="Sweden"]
wilcox.test(german_pc2, sweden_pc2)
#p=0.0309 W=5694

####How does load affect distinction
cor.test(g_s_load, pc1)

####Load vs PCs continued
load_pcoa <- ggplot(data=data.frame(bc_pcoa$vectors), aes(x=Axis.1, y=Axis.2)) + 
  geom_point(aes(color=as.numeric(as.character(log10(g_s_load))))) +
  scale_color_gradient2(name=expression(paste('log'[10],"(Tot. Load)", sep="")), low = "blue", mid = "white", high = "red") +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") + 
  theme(legend.key.width=unit(1,"cm"), axis.title.y = element_blank())
#ylab(paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))

pop_pcoa <-ggplot(data=data.frame(bc_pcoa$vectors), aes(x=Axis.1, y=Axis.2)) + 
  geom_point(aes(color=g_s_country)) +
  scale_color_discrete(name="Country") +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(paste(paste("Axis 1:", (100*round(perc_explained[1],3)), sep=" "), "%", sep="")) +
 theme(legend.key.width=unit(1,"cm"), axis.title.y = element_blank()) 
#  ylab(paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))

gA <- ggplotGrob(load_pcoa)
gB <- ggplotGrob(pop_pcoa)
gA$widths <- gB$widths

pdf(paste(paste("~/Dropbox/controlled_metagenomics/results_figures/pcoa_bray_curtis_", concat_name, sep=""), ".pdf", sep=""), family = "ArialMT")
grid.arrange(gA, gB, nrow = 2, left = paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))
dev.off()


########Jaccard
#PCoA on Jaccard
#jaccard <-as.matrix(vegdist(g_s_otus_filter, binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE, method="jaccard"))
#jaccard_pcoa=pcoa(jaccard, correction = "lingoes")
#bip=biplot(jaccard_pcoa, pch=20)
#perc_explained=jaccard_pcoa$values[,3]

plot_country <- ggplot(data=data.frame(jaccard_pcoa$vectors), aes(x=Axis.1, y=Axis.2)) + 
  geom_point(col= brewer.pal(n = 8, name = "Set1")[as.factor(g_s_country)]) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(paste(paste("Axis 1:", (100*round(perc_explained[1],3)), sep=" "), "%", sep="")) +
  ylab(paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))

#pdf("~/Dropbox/controlled_metagenomics/results_figures/pcoa_bray_curtis_g_s.pdf")
#plot_pcoa(bc_pcoa)
#dev.off()

####Do PC1 and PC2 distinguish between samples from germany vs sweden
pc1<-jaccard_pcoa$vectors[,1]
german_pc1<-pc1[g_s_country=="germany"]
sweden_pc1<-pc1[g_s_country=="Sweden"]
wilcox.test(german_pc1, sweden_pc1)
#p=0.0033 W=5367

pc2<-jaccard_pcoa$vectors[,2]
german_pc2<-pc2[g_s_country=="germany"]
sweden_pc2<-pc2[g_s_country=="Sweden"]
wilcox.test(german_pc2, sweden_pc2)
#p=0.0309 W=5694

####How does load affect distinction
#cor.test(g_s_load, pc1)

####Load vs PCs continued
load_pcoa <- ggplot(data=data.frame(bc_pcoa$vectors), aes(x=Axis.1, y=Axis.2)) + 
  geom_point(aes(color=as.numeric(as.character(log10(g_s_load))))) +
  scale_color_gradient2(name=expression(paste('log'[10],"(Tot. Load)", sep="")), low = "blue", mid = "white", high = "red") +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") + 
  theme(legend.key.width=unit(1,"cm"), axis.title.y = element_blank())
#ylab(paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))

pop_pcoa <-ggplot(data=data.frame(bc_pcoa$vectors), aes(x=Axis.1, y=Axis.2)) + 
  geom_point(aes(color=g_s_country)) +
  scale_color_discrete(name="Country") +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(paste(paste("Axis 1:", (100*round(perc_explained[1],3)), sep=" "), "%", sep="")) +
  theme(legend.key.width=unit(1,"cm"), axis.title.y = element_blank()) 
#  ylab(paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))

gA <- ggplotGrob(load_pcoa)
gB <- ggplotGrob(pop_pcoa)
gA$widths <- gB$widths

pdf(paste(paste("~/Dropbox/controlled_metagenomics/results_figures/pcoa_bray_curtis_g_s", concat_name, sep=""),".pdf",sep=""), family = "ArialMT")
grid.arrange(gA, gB, nrow = 2, left = paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))
dev.off()

##How is load associated with Shannon Diversity
shannon_div <- diversity(g_s_otus_filter)
simpson_div <- diversity(g_s_otus_filter, index="simpson")
J <- shannon_div/log(specnumber(g_s_otus_filter))
max_tax=apply(g_s_otus_filter, 1, max)

div_tog=bind_cols(load=g_s_load, shannon=shannon_div, simpson=simpson_div, J=J, max_tax=max_tax)

J_plot<-ggplot(data=div_tog, aes(x=load, y=J)) + 
  geom_point(pch=20) + 
  geom_smooth(method='lm', color='red', level=0) +
  theme_bw() + 
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Load (Total Microb. Cov./Plant Cov.)") +
  ylab("Evenness (Pielou's J)")

Shannon_div_plot<-ggplot(data=div_tog, aes(x=load, y=shannon)) + 
  geom_point(pch=20) + 
  geom_smooth(method='lm', color='red', level=0) +
  theme_bw() + 
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Load (Total Microb. Cov./Plant Cov.)") +
  ylab("Shannon Index (H')")

max_plot<-ggplot(data=div_tog, aes(x=load, y=max_tax*100)) + 
  geom_point(pch=20) + 
  geom_smooth(method='lm', color='red', level=0) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Load (Total Microb. Cov./Plant Cov.)") +
  ylab("Maximum Abundance (% of total)")

pdf(paste(paste("~/Dropbox/controlled_metagenomics/results_figures/div_load", concat_name, sep=""),".pdf",sep=""), family = "ArialMT")
plot_grid(Shannon_div_plot, max_plot, ncol=1)
dev.off()

##does load affect the overall evenness of the sample? Apparently not...?
cor.test(g_s_load, J)
#p-value = 0.6918, r=0.02604576
J_plot<-ggplot(J, g_s_load, pch=20)
cor.test(g_s_load, shannon_div)

##does load correlate with abundance of top taxonomical group

plot( log10(g_s_load), max_tax,pch=20)


#####BRAY CURTIS: Do these populations differ in who is there?


######Sweden only 
sweden_g_s <- g_s[which(g_s$country=="Sweden"),]
pc<-prcomp((sweden_filter), center = TRUE, scale. = TRUE, retx = TRUE)
ev <- pc$sdev^2

#Jari Oksanen explanation https://stat.ethz.ch/pipermail/r-help/2005-August/076610.html
prop_ev <- ev/sum(ev)
colfunc <- colorRampPalette(c("yellow", "purple"))
plot(pc$x[, 1], pc$x[, 2], col = colfunc(2)[as.factor(g_s_country)], main = "PCA", xlab = "PC1", ylab = "PC2", pch=20)

#x	if retx is true the value of the rotated data (the centred (and scaled if requested) data multiplied by the rotation matrix) is returned. Hence, cov(x) is the diagonal matrix diag(sdev^2). For the formula method, napredict() is applied to handle the treatment of values omitted by the na.action.

pc1 <- pc$x[,1]
pc2 <- pc$x[,2]
pc3 <- pc$x[,3]

}


#summary(lmer(pc1~(g_s_load)+(1|g_s_pop)))
#summary(lm(pc2~g_s_load+g_s_country))
#summary(lm(pc3~g_s_load+g_s_country))

#pca_vectors = data.frame(pc$x)
ggplot(data=pca_vectors, aes(x=PC1, y=PC2)) + 
  geom_point(col= brewer.pal(n = 8, name = "Set1")[as.factor(g_s_country)]) +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(paste(paste("Axis 1:", (100*round(prop_ev[1],3)), sep=" "), "%", sep="")) +
  ylab(paste(paste("Axis 2:", (100*round(prop_ev[2],3)), sep=" "), "%", sep=""))

ggplot(data=pca_vectors, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color=as.numeric(as.character(log10(g_s_load))))) +
  scale_color_gradient2(name=expression(paste('log'[10],"(Tot. Load)", sep="")),  low = "blue", mid = "white",high = "red") +
  theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") + 
  theme(legend.key.width=unit(1,"cm"), axis.title.y = element_blank())
#ylab(paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))
