### plotting georeferenced specimen collections by density

setwd("/Users/ethanlinck/Dropbox/miscellaneous_R/collection_priorities/")
library(ggplot2);library(ggmap);library(RColorBrewer);library(mapdata);library(maps)

localities <- read.csv("PNG_specimen_records.csv")
map <- map_data("worldHires", xlim=c(145,154), ylim=c(-10,0))
ggplot() + theme_bw() + 
  geom_path(data=map, aes(x=long, y=lat, group=group)) + coord_map() +
  stat_bin2d(data=localities, aes(x=decimallongitude, y=decimallatitude,fill = ..count..), binwidth = c(0.5,0.5)) +
  scale_fill_gradientn(colors = terrain.colors(10))
