##Bacteria
First we will look at those reads classified as bacterial. We will perform principal coordinate analysis on the bray-curtis disimilarity to assess how load relates to similarity of microbiome content.

```{r, echo=FALSE}

concat_name="_bacteria.csv"

g_s=read.csv(paste("~/Dropbox/controlled_metagenomics/results_figures/sweden_germany_combined", concat_name, sep=''), header = TRUE, row.names = 1, stringsAsFactors = FALSE)

g_s_tot = process_gs(g_s, subset=FALSE)

g_s_otus_filter = g_s_tot[[1]]

g_s_country = g_s_tot[[2]]

g_s_pop = g_s_tot[[3]]

g_s_load = g_s_tot[[4]]

bray_curtis <-as.matrix(vegdist(g_s_otus_filter, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE))
```

###Which microbe has the highest mean load across samples in a region?
```{r}
find_top_abundance(g_s_otus_filter)
```

###PCoA on bray-curtis
```{r}
bc_pcoa=pcoa(bray_curtis, correction = "lingoes")
perc_explained=bc_pcoa$values[,3]
```

####Load vs PCs continued
```{r}
load_pcoa <- ggplot(data=data.frame(bc_pcoa$vectors), aes(x=Axis.1, y=Axis.2)) + 
  geom_point(aes(color=log10(g_s_load)), cex = 2) +
  scale_color_gradient2(name = "", low = "blue", mid = "white", high = "red") +
  theme_light() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
  xlab("") + 
  theme(legend.key.width=unit(1,"cm"), axis.title.y = element_blank(), legend.text.align = 1)
# ylab(paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))

pop_pcoa <-ggplot(data=data.frame(bc_pcoa$vectors), aes(x=Axis.1, y=Axis.2)) + 
  geom_point(aes(color=g_s_country), cex = 2) +
  scale_color_discrete(name="Country") +
  theme_light() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
  xlab(paste(paste("Axis 1:", (100*round(perc_explained[1],3)), sep=" "), "%", sep="")) +
  theme(legend.key.width=unit(1,"cm"), axis.title.y = element_blank()) 
#  ylab(paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))

gA <- ggplotGrob(load_pcoa)
gB <- ggplotGrob(pop_pcoa)
gA$widths <- gB$widths

grid.arrange(gA, gB, nrow = 2, left = paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))

pdf(paste(paste("~/Dropbox/controlled_metagenomics/results_figures/pcoa_bray_curtis_", concat_name, sep=""), ".pdf", sep=""), family = "ArialMT", useDingbats=FALSE)
grid.arrange(gA, gB, nrow = 2, left = paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""), widths = 10)
dev.off()
```


###Do PC1 and PC2 distinguish between samples from germany vs sweden?
First PC1:
  ```{r, echo = FALSE}
pc1<-bc_pcoa$vectors[,1]
german_pc1<-pc1[g_s_country=="germany"]
sweden_pc1<-pc1[g_s_country=="Sweden"]
wilcox.test(german_pc1, sweden_pc1)
#p=0.1486 W=5367
```

###Next PC2:
```{r, echo = FALSE}
pc2<-bc_pcoa$vectors[,2]
german_pc2<-pc2[g_s_country=="germany"]
sweden_pc2<-pc2[g_s_country=="Sweden"]
wilcox.test(german_pc2, sweden_pc2)
#p=2.2e-16 W=8645
```

###What is the relationship between load and PC1?
```{r}
cor.test(g_s_load, pc1)
```


#PERMANOVA WITH BRAY-CURTIS DISTANCES
Here is some advice on nesting:
  https://stats.stackexchange.com/questions/188519/adonis-in-vegan-order-of-variables-or-use-of-strata

```{r}
#exclude NA
bc_comp = bray_curtis[which(!is.na(g_s_pop)), which(!is.na(g_s_pop))]
g_s_load_com = g_s_load[which(!is.na(g_s_pop))]
g_s_pop_com = g_s_pop[which(!is.na(g_s_pop))]
g_s_country_com = g_s_country[which(!is.na(g_s_pop))]
pop_permanova = adonis(bc_comp ~ g_s_pop_com / g_s_country_com + g_s_load_com, strata = g_s_pop_com )
country_permanova = adonis(bc_comp ~ g_s_load_com + g_s_country_com)
```







###How is load associated with Shannon Diversity

shannon_div <- vegan::diversity(g_s_otus_filter, index = "shannon")
simpson_div <- vegan::diversity(g_s_otus_filter, index="simpson")
J <- shannon_div/log(specnumber(g_s_otus_filter))
max_tax=apply(g_s_otus_filter, 1, max)
max_tax_names = colnames(g_s_otus_filter)[apply(g_s_otus_filter, 1, which.max)]

div_tog=bind_cols(load=g_s_load, shannon=shannon_div, simpson=simpson_div, J=J, max_tax=max_tax)

J_plot<-ggplot(data=div_tog, aes(x=log10(load), y=J)) + 
  geom_point(pch=20) + 
  geom_smooth(method = "glm", method.args = list(family = quasibinomial(link = 'logit')), color='red', level=0) +
  theme_light() + 
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression(paste('Load (log'[10],"(Total Microb. Cov./Plant Cov.))", sep=""))) +
  ylab("Evenness (Pielou's J)")

Shannon_div_plot<-ggplot(data=div_tog, aes(x=lload, y=shannon)) + 
  geom_point(pch=20) + 
  geom_smooth(method = "glm", method.args = list(family = binomial(link = 'logit')), color='red', level=0) +
  theme_light() + 
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression(paste(' Load (log'[10],"(Total Microb. Cov./Plant Cov.))", sep=""))) +
  ylab("Shannon Index (H')")

max_plot<-ggplot(data=div_tog, aes(x=log10(load), y=max_tax)) + 
  geom_point(pch=20) + 
  geom_smooth(method = "glm", method.args = list(family = quasibinomial(link = 'logit')), color='red', level=0) +
  theme_light() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression(paste(' Load (log'[10],"(Total Microb. Cov./Plant Cov.))", sep=""))) +
  ylab("Maximum Abundance (% of total)") + 
  scale_y_continuous(labels = function(x) x * 100) +
  geom_text(aes(label=ifelse(load>5,as.character(max_tax_names),'')),position=position_jitter(width=.1,height=.1))


plot_grid(J_plot, max_plot, ncol=1)

pdf(paste(paste("~/Dropbox/controlled_metagenomics/results_figures/div_load", concat_name, sep=""),".pdf",sep=""), family = "ArialMT")

plot_grid(J_plot, max_plot, ncol=1)

dev.off()
```

###Does load affect the overall evenness of the sample?
#```{r}
cor.test(g_s_load, J)

cor.test(g_s_load, shannon_div)
```

##Oomycetes

```{r, echo=FALSE}
concat_name="_oomycete.csv"

g_s=read.csv(paste("~/Dropbox/controlled_metagenomics/results_figures/sweden_germany_combined", concat_name, sep=''), header = TRUE, row.names = 1, stringsAsFactors = FALSE)

g_s_tot = process_gs(g_s, subset=FALSE)

g_s_otus_filter = g_s_tot[[1]]

g_s_country = g_s_tot[[2]]

g_s_pop = g_s_tot[[3]]

g_s_load = g_s_tot[[4]]

bray_curtis <-as.matrix(vegdist(g_s_otus_filter, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE))

#PCoA on bray-curtis
bc_pcoa=pcoa(bray_curtis, correction = "lingoes")

perc_explained=bc_pcoa$values[,3]

####Load vs PCs continued
load_pcoa <- ggplot(data=data.frame(bc_pcoa$vectors), aes(x=Axis.1, y=Axis.2)) + 
  geom_point(aes(color=log10(g_s_load)), cex = 2) +
  scale_color_gradient2(name = "", low = "blue", mid = "white", high = "red") +
  theme_light() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
  xlab("") + 
  theme(legend.key.width=unit(1,"cm"), axis.title.y = element_blank(), legend.text.align = 1)
# ylab(paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))

pop_pcoa <-ggplot(data=data.frame(bc_pcoa$vectors), aes(x=Axis.1, y=Axis.2)) + 
  geom_point(aes(color=g_s_country), cex = 2) +
  scale_color_discrete(name="Country") +
  theme_light() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 20)) +
  xlab(paste(paste("Axis 1:", (100*round(perc_explained[1],3)), sep=" "), "%", sep="")) +
  theme(legend.key.width=unit(1,"cm"), axis.title.y = element_blank()) 
#  ylab(paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))

gA <- ggplotGrob(load_pcoa)
gB <- ggplotGrob(pop_pcoa)
gA$widths <- gB$widths
gA <- ggplotGrob(load_pcoa)
gB <- ggplotGrob(pop_pcoa)
gA$widths <- gB$widths

grid.arrange(gA, gB, nrow = 2, left = paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))

pdf(paste(paste("~/Dropbox/controlled_metagenomics/results_figures/pcoa_bray_curtis_", concat_name, sep=""), ".pdf", sep=""), family = "ArialMT", useDingbats = FALSE)
grid.arrange(gA, gB, nrow = 2, left = paste(paste("Axis 2:", (100*round(perc_explained[2],3)), sep=" "), "%", sep=""))
dev.off()

```

###Which microbe has the highest load across samples?
```{r}
find_top_abundance(g_s_otus_filter)

```

###Do PC1 and PC2 distinguish between samples from germany vs sweden?
First PC1:
  ```{r, echo = FALSE}
pc1<-bc_pcoa$vectors[,1]
german_pc1<-pc1[g_s_country=="germany"]
sweden_pc1<-pc1[g_s_country=="Sweden"]
wilcox.test(german_pc1, sweden_pc1)
#p=0.0033 W=5367
```

###Next PC2:
```{r, echo = FALSE}
pc2<-bc_pcoa$vectors[,2]
german_pc2<-pc2[g_s_country=="germany"]
sweden_pc2<-pc2[g_s_country=="Sweden"]
wilcox.test(german_pc2, sweden_pc2)
#p=0.0309 W=5694
```

###What is the relationship between load and PC1?
```{r, echo = FALSE}
cor.test(g_s_load, pc1)
```


#PERMANOVA WITH BRAY-CURTIS DISTANCES
Here is some advice on nesting:
  https://stats.stackexchange.com/questions/188519/adonis-in-vegan-order-of-variables-or-use-of-strata

```{r}
#exclude NA
bc_comp = bray_curtis[which(!is.na(g_s_pop)), which(!is.na(g_s_pop))]
g_s_load_com = g_s_load[which(!is.na(g_s_pop))]
g_s_pop_com = g_s_pop[which(!is.na(g_s_pop))]
g_s_country_com = g_s_country[which(!is.na(g_s_pop))]
pop_permanova = adonis(bc_comp ~ g_s_pop_com / g_s_country_com + g_s_load_com, strata = g_s_pop_com )
country_permanova = adonis(bc_comp ~ g_s_load_com + g_s_country_com)
```























###How is load associated with Shannon Diversity
#```{r, echo=FALSE}
shannon_div <- vegan::diversity(g_s_otus_filter, index = "shannon")
simpson_div <- vegan::diversity(g_s_otus_filter, index="simpson")
J <- shannon_div/log(specnumber(g_s_otus_filter))
max_tax=apply(g_s_otus_filter, 1, max)
max_tax_names = colnames(g_s_otus_filter)[apply(g_s_otus_filter, 1, which.max)]
div_tog=bind_cols(load=g_s_load, shannon=shannon_div, simpson=simpson_div, J=J, max_tax=max_tax)

J_plot<-ggplot(data=div_tog, aes(x=log10(load), y=J)) + 
  geom_point(pch=20) + 
  geom_smooth(method = "glm", 
              method.args = list(family = quasibinomial(link = 'logit')), color='red', level=0) +
  theme_light() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Load (Total Microb. Cov./Plant Cov.)") +
  ylab("Evenness (Pielou's J)")

Shannon_div_plot<-ggplot(data=div_tog, aes(x=log10(load), y=shannon)) + 
  geom_point(pch=20) + 
  geom_smooth(method='glm', method.args = list(family = quasibinomial(link = 'logit')), color='red', level=0) +
  theme_light() + 
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression(paste(' Load (log'[10],"(Total Microb. Cov./Plant Cov.))", sep=""))) +
  ylab("Shannon Index (H')")

max_plot<-ggplot(data=div_tog, aes(x=log10(load), y=max_tax)) + 
  geom_point(pch=20) + 
  geom_smooth(method = "glm", 
              method.args = list(family = quasibinomial(link = 'logit')), color='red', level=0) +
  theme_light() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression(paste(' Load (log'[10],"(Total Microb. Cov./Plant Cov.))", sep=""))) +
  ylab("Maximum Abundance (% of total)")  + 
  scale_y_continuous(labels = function(x) x * 100) +
  geom_text(aes(label=ifelse(load>0.05,as.character(max_tax_names),'')),position=position_jitter(width=.3,height=.3))

plot_grid(J_plot, max_plot, ncol=1)

pdf(paste(paste("~/Dropbox/controlled_metagenomics/results_figures/div_load", concat_name, sep=""),".pdf",sep=""), family = "ArialMT")
plot_grid(J_plot, max_plot, ncol=1)
dev.off()
```

###Does load affect the overall evenness of the sample?
#```{r, echo = FALSE}
cor.test(g_s_load, J)

cor.test(g_s_load, shannon_div)
```
