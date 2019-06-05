#the goal of this script is to get the metagenome tables into a shape where they can be analyzed downstream. This is used in the publication
#bacteria sweden
#concat_name="bacteria.csv"
for(concat_name in c("bacteria.csv", "fungi.csv", "oomycete.csv", "virus.csv")){
  process_class(concat_name)
}

process_class<-function(concat_name){
  #bacteria germany ***NEED TO FIX DUPLICATES***
  german_samples=read.table("~/Dropbox/germany_pathogen_collections/sample_data/plate_sample_locations/sample_infoFinal_2018.txt", sep="\t", header=T, stringsAsFactors = FALSE)
  keep=german_samples #german_samples[-which(duplicated(german_samples$uniqueID)),]
  meta_sweden=read.table(paste("~/work_main/abt6_projects9/metagenomic_controlled/data/processed_reads/swedish_samples/meta_family_corrected_per_plant_", concat_name, sep=""), sep=",",header=T, row.names = 1, stringsAsFactors = FALSE)
  meta_german1=read.table(paste("~/work_main/abt6_projects9/metagenomic_controlled/data/processed_reads/german_samples/meta_family_corrected_per_plant_", concat_name, sep=""), sep=",",header=T, row.names = 1, stringsAsFactors = FALSE)

  colnames(meta_sweden)=gsub( "_RunId0094_LaneId7","", colnames(meta_sweden))
  population=lapply(strsplit(colnames(meta_sweden), "\\."), "[[",1)
  sample=unlist(strsplit(colnames(meta_sweden), '\\.'))[c(FALSE, TRUE)]
  meta_sweden["load",]=colSums(meta_sweden, na.rm=T)
  meta_sweden["population",]=population
  meta_sweden["country",]="Sweden"

  #remove controls
  meta_sweden=meta_sweden[,colnames(meta_sweden)!=c("Control_1")]
  meta_sweden=meta_sweden[,colnames(meta_sweden)!=c("Control_2")]
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
  g_s[is.na(g_s)]<-0

  write.csv(g_s, paste("~/Dropbox/controlled_metagenomics/results_figures/sweden_germany_combined_", concat_name, sep=""), row.names = TRUE, quote = FALSE)

  write.csv(t_sweden, paste("~/Dropbox/controlled_metagenomics/results_figures/sweden_only_", concat_name, sep=""),  row.names = TRUE, quote = FALSE)

  write.csv(t_german, paste("~/Dropbox/controlled_metagenomics/results_figures/germany_only_", concat_name, sep=""),  row.names = TRUE, quote = FALSE)
}


