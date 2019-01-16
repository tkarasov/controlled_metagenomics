library(ggmap)
scalebar = function(x,y,w,n,d, units="km"){
  # x,y = lower left coordinate of bar
  # w = width of bar
  # n = number of divisions on bar
  # d = distance along each division
  
  bar = data.frame( 
    xmin = seq(0.0, n*d, by=d) + x,
    xmax = seq(0.0, n*d, by=d) + x + d,
    ymin = y,
    ymax = y+w,
    z = rep(c(1,0),n)[1:(n+1)],
    fill.col = rep(c("black","white"),n)[1:(n+1)])
  
  labs = data.frame(
    xlab = c(seq(0.0, (n+1)*d, by=d) + x, x), 
    ylab = c(rep(y+w*1.5, n+2), y+3*w),
    text = c(as.character(seq(0.0, ((n+1)*d)/1000, by=d/1000)), units)
  )
  list(bar, labs)
}

#germany
locations=read.table("~/Dropbox/germany_pathogen_collections/data_files_rmarkdown/location_sampling_lat_long.txt",header=T)
#sb = scalebar(9.05,48.45, 0.005, , 0.0143, "km" )

pdf("~/Dropbox/germany_pathogen_collections/data_files_rmarkdown/locations_collections.pdf")
qmplot(Longitude, Latitude, data = locations, extent = "panel")+geom_point(color="RED", cex=4)+geom_rect(data=sb[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F, show.legend = F,  color = "black", fill = sb[[1]]$fill.col)+geom_text(data=sb[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F) 
dev.off()

#sweden
locations=read.table("~/Dropbox/controlled_metagenomics/Swedish_samples/swedish_sites.txt",header=T)
sb = scalebar(9.05,48.45, 0.005,5, 0.0143, "km" )

pdf("~/Dropbox/controlled_metagenomics/results_figures/swedish_map.pdf")
qmplot(Longitude, Latitude, data = locations, extent = "panel")+geom_point(color="RED", cex=4)+geom_rect(data=sb[[1]], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z), inherit.aes=F, show.legend = F,  color = "black", fill = sb[[1]]$fill.col)+geom_text(data=sb[[2]], aes(x=xlab, y=ylab, label=text), inherit.aes=F, show.legend = F) 
dev.off()
