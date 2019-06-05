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

plot_side_by_side<-function(out, colors){
  max_y_axis=max(max(out$germany[,'value']), max(out$sweden[,'value'])) + 0.2*max(max(out$germany[,'value']), max(out$sweden[,'value']))
  
  b_g = ggplot(data=out$germany, aes(x=Load, y=value, fill=Family))+ geom_bar(aes(), stat="identity", position="stack") +
    scale_fill_manual(values = colors) +
    theme(legend.position="bottom", panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank()) + 
    guides(fill=guide_legend(nrow=5)) +
    xlab("Plant Individuals") +
    theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_blank(), legend.title = element_blank())+
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_y_axis)) +
    annotate("text",  x=Inf, y = Inf, label = "Germany", vjust=7, hjust=2, cex = 5)
  
  s_g = ggplot(data=out$sweden, aes(x=Load, y=value, fill=Family))+ geom_bar(aes(), stat="identity", position="stack") +
    scale_fill_manual(values = colors) +
    theme(legend.position="bottom", legend.title = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank()) + 
    guides(fill = guide_legend(nrow = 5)) +
    xlab("Plant Individuals") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_y_axis)) +
    annotate("text",  x = Inf, y = Inf, label = "Sweden", vjust=7, hjust=2, cex = 5) 
  
  grobs <- ggplotGrob(s_g)$grobs
  legend<-grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  pgrid<-plot_grid(b_g+theme(legend.position="none" , panel.border = element_rect(colour = "black", fill=NA, size=1)), s_g + theme(legend.position="none") +
                     theme(axis.title.y=element_blank()), nrow = 2)
  p <- plot_grid(pgrid, legend, nrow = 2, ncol=1, rel_heights = c(2, .45))
  
  p
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
  #g_s_otus <-g_s[,which(apply((g_s), 1, var)!=0)]
  #first samples that are missing
  g_s_otus <-g_s_otus[which(apply((g_s_otus), 1, var)!=0),]
  #Next taxa that are absent
  g_s_otus <-g_s_otus[which(apply(t(g_s_otus), 2, var)!=0),]
  #remove families that are 
  keep_0.001 <- which(colSums(g_s_otus)/dim(g_s_otus)[1]>0.00001)
  g_s_otus_filter <- g_s_otus[,keep_0.001]
  #now make compositional
  g_s_otus_filter <- g_s_otus_filter/rowSums(g_s_otus)
  together<-list(g_s_otus_filter, g_s_country, g_s_pop, g_s_load)
  return(together)
}
