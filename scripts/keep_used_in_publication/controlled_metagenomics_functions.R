
theme_tufte_revised <- function(base_size = 11, base_family = "ArialMT", ticks = TRUE) {
  
  ret <- theme_bw(base_family = base_family, base_size = base_size) + 
    theme(
      axis.line = element_line(color = 'black'),
      axis.title.x = element_text(vjust = -0.3), 
      axis.title.y = element_text(vjust = 0.8),
      legend.background = element_blank(), 
      legend.key = element_blank(), 
      legend.title = element_text(face="plain"),
      panel.background = element_blank(), 
      panel.border = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_blank(),
      strip.background = element_blank()
    )
  
  if (!ticks) {
    ret <- ret + theme(axis.ticks = element_blank())
  }
  
  ret
} 

meta_table <- function(level){
  concat_name = paste(level,".csv",sep="")
  meta_both = read.table(paste("~/Dropbox/controlled_metagenomics/data/sweden_germany_combined_nodup_", concat_name, sep=""), sep=",",header=T, row.names = 1)
  meta_sweden0 = meta_both %>% filter(country == "Sweden")
  meta_german0 = meta_both %>% filter(country == "germany")
  meta_sweden = meta_sweden0 %>% select(-c("country", "population", "load"))
  meta_german = meta_german0 %>% select(-c("country", "population", "load"))
  top_swed = sort(colSums(meta_sweden, na.rm=TRUE), decreasing=TRUE )[1:5]
  top_germ = sort(colSums(meta_german, na.rm=TRUE), decreasing=TRUE )[1:5]
  together = unique(c(names(top_germ), names(top_swed)))
  
  germany_processed = process_microb_melt(meta_german, together)
  germany_processed$value = as.numeric(as.character(germany_processed$value))
  germany_processed$Family = as.factor(germany_processed$Family)
  
  sweden_processed = process_microb_melt(meta_sweden,together)
  sweden_processed$value = as.numeric(as.character(sweden_processed$value))
  sweden_processed$Family = as.character(sweden_processed$Family)
  out<-list()
  out$germany = germany_processed
  out$sweden = sweden_processed
  return(out)
}

process_microb_melt<-function(meta, together){
  meta_microbiome = meta[,together]
  rest = rowSums(meta[,colnames(meta) %ni% together ], na.rm=TRUE)
  meta_microbiome$Other = rest
  tot = rowSums(meta_microbiome, na.rm=TRUE)
  meta_microbiome = as.data.frame(t(meta_microbiome))
  meta_microbiome$Family = rownames(meta_microbiome)
  microb_melt = melt(meta_microbiome, id=c("Family"))
  microb_melt$Load = as.factor(-1*tot[microb_melt$variable])
  return(microb_melt)
}

plot_side_by_side<-function(out, colors){
  max_y_axis=max(max(out$germany[,'value']), max(out$sweden[,'value'])) + 0.2*max(max(out$germany[,'value']), max(out$sweden[,'value']))

  b_g = ggplot(data=out$germany, aes(x=Load, y=value, fill=Family))+ geom_bar(aes(), stat="identity", position="stack") +
    scale_fill_manual(values = colors) +
    theme(legend.position="bottom", 
          panel.background = element_blank(), 
          axis.text.x = element_blank(), 
          axis.title.x=element_blank() , 
          axis.ticks.x=element_blank(), 
          axis.title.y = element_blank(), 
          legend.title = element_blank(), 
          panel.border = element_rect(color = "gray50")) + 
          #base_size = 11, base_family = "",
          #base_line_size = base_size/22, base_rect_size = base_size/22) + 
    guides(fill=guide_legend(nrow=5)) +
    xlab("Plant Individuals") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_y_axis)) +
    annotate("text",  x=Inf, y = Inf, label = "Germany", vjust=7, hjust=2, cex = 5) 
  
  s_g = ggplot(data=out$sweden, aes(x=Load, y=value, fill=Family))+ geom_bar(aes(), stat="identity", position="stack") +
    scale_fill_manual(values = colors) +
    theme(legend.position="bottom", panel.background = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), axis.line = element_line(color = "grey10") , axis.ticks.x=element_blank(), axis.title.y = element_blank(), legend.title = element_blank()) + 
    guides(fill=guide_legend(nrow=3)) +
    xlab("Plant Individuals") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_y_axis)) +
    annotate("text",  x=Inf, y = Inf, label = "Sweden", vjust=7, hjust=2, cex = 5) 
  
  grobs <- ggplotGrob(s_g)$grobs
  
  legend<-grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  
  pgrid <- plot_grid(b_g + theme(legend.position="none", panel.border = element_rect(colour = "grey10", fill=NA, size=1)), s_g + theme(legend.position="none", panel.border = element_rect(colour = "grey0", fill=NA, size=1)) + theme(axis.title.y=element_blank()), nrow = 2)
  p <- plot_grid(pgrid, legend, nrow = 2, ncol=1, rel_heights = c(2, .45))
  
  p
}


reduce_melted <- function(out, cutoff){
  #Cutoff is bottom percentile
  german_cutoff = distinct(out$germany[,c("variable", "Load")])
  german_cutoff = german_cutoff[order(german_cutoff$Load, decreasing = T),]
  german_cutoff = german_cutoff[1:round(cutoff*dim(german_cutoff)[1]),]
  g_cut = out$germany[which(out$germany$variable %in% german_cutoff$variable),]
  
  sweden_cutoff = distinct(out$sweden[,c("variable", "Load")])
  sweden_cutoff = sweden_cutoff[order(sweden_cutoff$Load, decreasing = T),]
  sweden_cutoff = sweden_cutoff[1:round(cutoff*dim(sweden_cutoff)[1]),]
  s_cut = out$sweden[which(out$sweden$variable %in% sweden_cutoff$variable),]
  
  out$germany = g_cut
  out$sweden = s_cut
  
  return(out)
}


find_top_abundance <- function(g_s_otus_filter){
  #Which families have the highest loads in the different countries?
  swedish = g_s_otus_filter[which(g_s_country == "Sweden"),]
  swede_top = colnames(swedish)[apply(swedish, 1, which.max)]
  tot_top = max(table(swede_top))/sum(table(swede_top))
  ident_top = which(table(swede_top)==max(table(swede_top)))
  max_sweden = names(which.max(apply(swedish, 2, mean)))
  german = g_s_otus_filter[which(g_s_country == "germany"),]
  german_top = colnames(german)[apply(german, 1, which.max)]
  tot_top = max(table(german_top))/sum(table(german_top))
  ident_top = which(table(german_top)==max(table(german_top)))
  max_german = names(which.max(apply(german, 2, mean)))
  
  return(c(max_german, max_sweden))
}

process_gs<-function(g_s, subset=FALSE, classifier=cf, classifier_identity=cfi){
  #This function processes an original input file and returns the compositionally corrected table
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

  #first samples that are missing
  g_s_otus <-g_s_otus[which(apply((g_s_otus), 1, var)!=0),]
  
  #Next taxa that are absent
  g_s_otus <-g_s_otus[which(apply(t(g_s_otus), 2, var)!=0),]
  
  #remove families that are found at an average
  keep_0.001 <- which(colSums(g_s_otus)/dim(g_s_otus)[1]>0.00001)
  g_s_otus_filter <- g_s_otus[,keep_0.001]
  
  #now make compositional
  g_s_otus_filter_comp <- g_s_otus_filter/rowSums(g_s_otus)
  together<-list(g_s_otus_filter_comp, g_s_country, g_s_pop, g_s_load, g_s_otus_filter)
  return(together)
}

process_gs_keep_all <-function(g_s, subset=FALSE, classifier=cf, classifier_identity=cfi){
  #This function processes an original input file and returns the compositionally corrected table
  #function that takes gs
  if(subset==TRUE){
    g_temp<-g_s[which(g_s[[cf]]==cfi),]
  }
  else{
    g_temp<-g_s
  }
  g_s_otus <- g_temp[,which(colnames(g_temp)%ni%c("population", "country", "load"))]
  g_s_load <-g_temp$load
  g_s_pop <-g_temp$population
  g_s_country <-g_temp$country
  
  #first samples that are missing
  g_s_otus <-g_s_otus[which(apply((g_s_otus), 1, var)!=0),]
  
  #Next taxa that are absent
  g_s_otus <-g_s_otus[which(apply(t(g_s_otus), 2, var)!=0),]
  
  #remove families that are found at an average greater than 0.00001
  keep_0.001 <- which(colSums(g_s_otus)/dim(g_s_otus)[1]>0.00001)
  g_s_otus_filter <- g_s_otus[,keep_0.001]
  
  #now make compositional
  g_s_otus_filter <- g_s_otus_filter/rowSums(g_s_otus)
  together<-list(g_s_otus_filter, g_s_country, g_s_pop, g_s_load)
  return(together)
}

merge_phyloseq<-function(g_s){
  rownames(g_s) = g_s$Genotype
  g_s = g_s[,which(colnames(g_s)!="Genotype")]
  #g_s_otus <- g_s[,which(colnames(g_s)%ni%c("population", "country", "load", "Genotype"))]
  g_s_tot = process_gs(g_s, subset=FALSE)
  g_s_otus_comp = g_s_tot[[1]]
  g_s_country = g_s_tot[[2]]
  g_s_pop = g_s_tot[[3]]
  g_s_load = g_s_tot[[4]]
  
  #g_s_filter is the data that is not compositional but instead coverage per each genome
  g_s_filter = apply(round(g_s_tot[[5]] * 12000000), 2, as.integer)
  g_s_filter[is.na(g_s_filter)] = 0
  OTU = (otu_table(g_s_filter, taxa_are_rows = FALSE))
  row.names(OTU) = row.names(g_s_tot[[5]])
  
  #Prepare sample data
  sampledata = sample_data(data.frame(Population = g_s_pop, Country = g_s_country, load = g_s_load, row.names = row.names(OTU)))
  physeq = phyloseq(otu_table=OTU, sam_data=sampledata)
  return(physeq)
}

deseq_variance_stab<-function(physeq, deseq_data){
  trial.DEseq = estimateSizeFactors(deseq_data)
  trial.DEseq = estimateDispersions(trial.DEseq)
  trial.vst = getVarianceStabilizedData(trial.DEseq)
  dim(trial.vst)
  
  # Save the untransformed data as a separate variable so you can go back to it
  trial.phyloseq0 = physeq
  
  # add the varience stabilized otu numbers into the dataset:
  otu_table(trial.phyloseq0) <- otu_table(trial.vst, taxa_are_rows = TRUE)
  
  return(trial.phyloseq0)
}


