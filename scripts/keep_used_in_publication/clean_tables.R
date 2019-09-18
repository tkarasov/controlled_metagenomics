#the goal of this script is to get the metagenome tables into a shape where they can be analyzed downstream. This is used in the publication
#bacteria sweden
#concat_name="bacteria.csv"

process_class<-function(concat_name){
  #bacteria germany ***NEED TO FIX DUPLICATES***
  german_samples=read.table("~/Dropbox/germany_pathogen_collections/sample_data/plate_sample_locations/sample_infoFinal_2018.txt", sep="\t", header=T, stringsAsFactors = FALSE)
  keep=german_samples #german_samples[-which(duplicated(german_samples$uniqueID)),]
  
  meta_sweden=read.table(paste("~/Dropbox/controlled_metagenomics/data/swedish_meta_family_corrected_per_plant_", concat_name, sep=""), sep=",",header=T, row.names = 1, stringsAsFactors = FALSE)
  
  meta_german1=read.table(paste("~/Dropbox/controlled_metagenomics/data/german_meta_family_corrected_per_plant_", concat_name, sep=""), sep=",",header=T, row.names = 1, stringsAsFactors = FALSE)

  colnames(meta_sweden)=gsub( "_RunId0094_LaneId7","", colnames(meta_sweden))
  population=lapply(strsplit(colnames(meta_sweden), "\\."), "[[",1)
  sample=unlist(strsplit(colnames(meta_sweden), '\\.'))[c(FALSE, TRUE)]
  meta_sweden["load",]=colSums(meta_sweden, na.rm=T)
  meta_sweden["population",]=population
  meta_sweden["country",]="Sweden"

  #remove controls
  meta_sweden=meta_sweden[,colnames(meta_sweden)!=c("Control_1")]
  meta_sweden=meta_sweden[,colnames(meta_sweden)!=c("Control_2")]
  #remove two bad swedish samples
  meta_sweden = meta_sweden[,which(colnames(meta_sweden) %ni% c("Adal.5", "ULL.5"))]
  t_sweden=data.frame(t(meta_sweden))
  temp=t_sweden[,!colnames(t_sweden)%in%c("country", "population")]
  t_sweden[,!colnames(t_sweden)%in%c("country", "population")]=mutate_all(temp, function(x) as.numeric(as.character(x)))

  meta_german=meta_german1[,keep$metagenome_identifier]
  meta_german["load",]=colSums(meta_german, na.rm=T)
  meta_german["country",]="germany"
  which_keep=keep[keep$metagenome_identifier%in%colnames(meta_german),]
  meta_german=meta_german[-which(duplicated(which_keep$uniqueID)),]
  population=german_samples[match(german_samples$metagenome_identifier,colnames(meta_german)),]$site
  meta_german['population',]=population
  t_german=data.frame(t(meta_german), stringsAsFactors = FALSE)
  temp=t_german[,!colnames(t_german)%in%c("country", "population")]
  t_german[,!colnames(t_german)%in%c("country", "population")]=mutate_all(temp, function(x) as.numeric(as.character(x)))

  #german_swedish
  t_sweden$name=rownames(t_sweden)
  t_german$name=rownames(t_german)
  g_s=merge(t_sweden, t_german, all=TRUE)

  rownames(g_s)=g_s$name
  g_s<-g_s[,which(colnames(g_s)!="name")]
  
  #assign zero to NA. Factors stay as NA
  g_s[is.na(g_s)]<-0

  write.csv(g_s, paste("~/Dropbox/controlled_metagenomics/data/sweden_germany_combined_", concat_name, sep=""), row.names = TRUE, quote = FALSE)

  write.csv(t_sweden, paste("~/Dropbox/controlled_metagenomics/data/sweden_only_", concat_name, sep=""),  row.names = TRUE, quote = FALSE)

  write.csv(t_german, paste("~/Dropbox/controlled_metagenomics/data/germany_only_", concat_name, sep=""),  row.names = TRUE, quote = FALSE)
}

process_untransformed_class<-function(concat_name){
  #This process the original output
  #bacteria germany ***NEED TO FIX DUPLICATES*** DONE
  german_samples=read.table("~/Dropbox/germany_pathogen_collections/sample_data/plate_sample_locations/sample_infoFinal_2018.txt", sep="\t", header=T, stringsAsFactors = FALSE)
  keep=german_samples #german_samples[-which(duplicated(german_samples$uniqueID)),]
  
  meta_sweden=read.table(paste("~/Dropbox/controlled_metagenomics/data/swedish_centrifuge_metagenome_table_", concat_name, sep=""), sep="\t",header=T, row.names = 1, stringsAsFactors = FALSE)
  
  meta_german1=read.table(paste("~/Dropbox/controlled_metagenomics/data/german_centrifuge_metagenome_table_", concat_name, sep=""), sep="\t",header=T, row.names = 1, stringsAsFactors = FALSE)
  
  colnames(meta_sweden)=gsub( "_RunId0094_LaneId7","", colnames(meta_sweden))
  population=lapply(strsplit(colnames(meta_sweden), "\\."), "[[",1)
  sample=unlist(strsplit(colnames(meta_sweden), '\\.'))[c(FALSE, TRUE)]
  meta_sweden["load",]=colSums(meta_sweden, na.rm=T)
  meta_sweden["population",]=population
  meta_sweden["country",]="Sweden"
  
  #remove controls
  meta_sweden=meta_sweden[,colnames(meta_sweden)!=c("Control_1")]
  meta_sweden=meta_sweden[,colnames(meta_sweden)!=c("Control_2")]
  colnames(meta_sweden) = gsub(".R1.fq.report", "",gsub("X.ebio.abt6_projects9.metagenomic_controlled.data.processed_reads.swedish_samples.centrifuge_output.", "", colnames(meta_sweden)))
  t_sweden=data.frame(t(meta_sweden))
  temp=t_sweden[,!colnames(t_sweden)%in%c("country", "population")]
  t_sweden[,!colnames(t_sweden)%in%c("country", "population")]=mutate_all(temp, function(x) as.numeric(as.character(x)))
  
  colnames(meta_german1) = gsub(".R1.fq.report", "",gsub("X.ebio.abt6_projects9.metagenomic_controlled.data.processed_reads.german_samples.centrifuge_output.", "", colnames(meta_german1)))
  meta_german=meta_german1[,keep$metagenome_identifier]
  meta_german["load",]=colSums(meta_german, na.rm=T)
  meta_german["country",]="germany"
  which_keep=keep[keep$metagenome_identifier%in%colnames(meta_german),]
  meta_german=meta_german[-which(duplicated(which_keep$uniqueID)),]
  population=german_samples[match(german_samples$metagenome_identifier,colnames(meta_german)),]$site
  meta_german['population',]=population
  t_german=data.frame(t(meta_german), stringsAsFactors = FALSE)
  temp=t_german[,!colnames(t_german)%in%c("country", "population")]
  t_german[,!colnames(t_german)%in%c("country", "population")]=mutate_all(temp, function(x) as.numeric(as.character(x)))
  
  #german_swedish
  t_sweden$name=rownames(t_sweden)
  t_german$name=rownames(t_german)
  g_s=merge(t_sweden, t_german, all=TRUE)
  
  rownames(g_s)=g_s$name
  g_s<-g_s[,which(colnames(g_s)!="name")]
  
  #assign zero to NA. Factors stay as NA
  g_s[is.na(g_s)]<-0
  
  write.csv(g_s, paste("~/Dropbox/controlled_metagenomics/data/original_sweden_germany_combined_", concat_name, sep=""), row.names = TRUE, quote = FALSE)
  
  write.csv(t_sweden, paste("~/Dropbox/controlled_metagenomics/data/original_sweden_only_", concat_name, sep=""),  row.names = TRUE, quote = FALSE)
  
  write.csv(t_german, paste("~/Dropbox/controlled_metagenomics/data/original_germany_only_", concat_name, sep=""),  row.names = TRUE, quote = FALSE)
}

for(concat_name in c("bacteria.csv", "fungi.csv", "oomycete.csv")){
  process_class(concat_name)
}

for(concat_name in c("bac.txt", "fungi.txt", "oom.txt")){
  process_untransformed_class(concat_name)
}

german_samples = read.table("~/Dropbox/germany_pathogen_collections/sample_data/plate_sample_locations/sample_infoFinal_2018.txt", sep="\t", header=T)

#duplicated samples (samples collected or sequenced twice)
keep = german_samples[-which(duplicated(german_samples$uniqueID)),]
no_keep = as.character(german_samples[which(duplicated(german_samples$uniqueID)),]$metagenome_identifier)
#These have too few reads: NextMet175, NextMet193, NextMet194, NextMet55
#These are the controls
#duplicated: NextMet24, NextMet50, NextMet66, NextMet70, NextMet83, NextMet85, NextMet87, NextMet96, NextMet119, NextMet131, NextMet142, NextMet144

keep = keep[which(keep$uniqueID %ni% c("NextMet175", "NextMet193", "NextMet194", "NextMet55")),]
no_keep = c(no_keep, c("NextMet175", "NextMet193", "NextMet194", "NextMet55")) 

#Let's get rid of bad columns for bacteria, oomycetes and fungi for the transformed and untransformed

remove_bad <-function(concat_name, original = FALSE){
  if(original == FALSE){
  g_s_all = read.csv(paste("~/Dropbox/controlled_metagenomics/data/sweden_germany_combined", concat_name, sep=''), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  g_s_keep = g_s_all[rownames(g_s_all) %ni% no_keep,]
  write.csv(g_s_keep, paste("~/Dropbox/controlled_metagenomics/data/sweden_germany_combined_nodup", concat_name, sep=''),row.names = TRUE)
  }
  if(original == TRUE){
    g_s_all = read.csv(paste("~/Dropbox/controlled_metagenomics/data/original_sweden_germany_combined", concat_name, sep=''), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    g_s_keep = g_s_all[rownames(g_s_all) %ni% no_keep,]
    write.csv(g_s_keep, paste("~/Dropbox/controlled_metagenomics/data/original_sweden_germany_combined_nodup", concat_name, sep=''),row.names = TRUE)
  }
  }

for(concat_name in c("_bacteria.csv", "_fungi.csv", "_oomycete.csv", "_virus.csv")){
  remove_bad(concat_name, original = FALSE)
}

for(concat_name in c("_bac.txt", "_oom.txt")){
  remove_bad(concat_name, original = TRUE)
}
