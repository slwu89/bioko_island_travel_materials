####
# Clean Raw Data - PR map draws
#
# Load in 100 map draws as well as mean + sd
# Format cleaned output as areaID, mean.PR, sd.PR, draw.1, draw.2, etc
#
# August 12, 2019
#
####

library(data.table)
library(sp)
library(raster)
library(rgdal)
library(here)
library(ggplot2)
here()


centroids <- rgdal::readOGR(here("data/raw/mapArea_centroids/mapareagrid_centroids.shp"))
dat <- raster::raster("data/raw/Output_popMCD/Prevalence_posterior_realization_1.tif")
# do not attempt to reconcile the two by setting the crs to be the same
# somehow, because the extents are so different, this really does mess up a lot
# but if I don't do anything before calling the extract

# Use the focal() function to deal with the edge effects
# https://gis.stackexchange.com/questions/187410/how-do-i-get-rid-of-edge-effects-while-using-focal-in-r-to-smooth-a-raster
dat.1 <- focal(dat, w = matrix(1, nc = 3, nr = 3), fun = mean, na.rm = TRUE)

pfpr.dt = data.table(
  coordinates(centroids), # lat-long coordinates of the centroids
  areaId = centroids$FID, # this is the same as areaId
  pfpr = extract(dat, centroids), # extract out the values of the raster map at the centroids of each pixel
  pfpr.edges.1 = extract(dat.1, centroids) # extract out values of the edge-corrected raster map at centroids
)
# unfortunately, we are missing some data around the edges when we apply the extraction
# and try to find the value of the raster at all of the centroids;
# instead, we build out the edges using the focal function() and sample from there
# not elegant, but here's how we fill in the missing edges on the same surface as before
pfpr.dt[is.na(pfpr) & !is.na(pfpr.edges.1)]$pfpr <- pfpr.dt[is.na(pfpr) & !is.na(pfpr.edges.1)]$pfpr.edges.1
pfpr.dt$pfpr.edges.1 <- NULL



bioko <- rgdal::readOGR("/Users/dtcitron/Documents/MASH/MASH-Main/MASH-dev/DanielCitron/Bioko_Island_Cluster_Simulations/BI_maps/bioko", "bioko_admin1")
areas_inh <- rgdal::readOGR("/Users/dtcitron/Documents/MASH/MASH-Main/MASH-dev/DanielCitron/Bioko_Island_Cluster_Simulations/BI_maps/areas_inh", "areas_inh")
areasf <- fortify(areas_inh, region = "areaId")

# and here's how we make a map
area.data = merge(areasf, pfpr.dt, by.x = "id", by.y = "areaId", all=FALSE)
plot.data <- area.data[order(area.data$order), ]
p1 = ggplot(data = plot.data, aes(x=long, y=lat, group = group))
p2 = p1 + geom_polygon(data = bioko, aes(x = long, y = lat, group = group), color = "black", fill="grey", size = 0.25)
map <- p2 + geom_polygon(data = plot.data, aes(x = long, y = lat, group = group, fill = pfpr), color = NA, size = 0.25) +
 scale_fill_gradient(name="pfpr", low="yellow", high="red", limits=c(0,.6)) +
 geom_polygon(data = bioko, aes(x = long, y = lat, group = group), color = "black", fill=NA, size = 0.25) +
 theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(),
       axis.title.x=element_blank(),
       axis.title.y=element_blank(), panel.background=element_blank(), legend.position=c(0.2, 0.8))
map

# Construct cleaned PfPR map data set ####
# Now the trick will be to read in all 100 files and store the extracted values for pfpr
# Everything will be indexed by areaId

# create blank data table
dt <- data.table(areaId = integer(),
                 pfpr = numeric(),
                 draw = integer()
)

for(j in 1:100){
  filename <- paste0(here("data/raw/Output_popMCD/"), "Prevalence_posterior_realization_",as.numeric(j),".tif")
  dat <- raster::raster(filename)
  dat.1 <- focal(dat, w = matrix(1, nc = 5, nr = 5), fun = mean, na.rm = TRUE)
  
  pfpr.dt = data.table(
    #coordinates(centroids), # lat-long coordinates of the centroids
    areaId = centroids$FID, # this is the same as areaId
    pfpr = extract(dat, centroids), # extract out the values of the raster map at the centroids of each pixel
    pfpr.edges.1 = extract(dat.1, centroids) # extract out values of the edge-corrected raster map at centroids
  )
  # unfortunately, we are missing some data around the edges when we apply the extraction
  # and try to find the value of the raster at all of the centroids;
  # instead, we build out the edges using the focal function() and sample from there
  # not elegant, but here's how we fill in the missing edges on the same surface as before
  pfpr.dt[is.na(pfpr) & !is.na(pfpr.edges.1)]$pfpr <- pfpr.dt[is.na(pfpr) & !is.na(pfpr.edges.1)]$pfpr.edges.1
  pfpr.dt$pfpr.edges.1 <- NULL
  pfpr.dt$draw <- j
  
  dt <- rbind(dt, pfpr.dt)
}
dt[, draw := paste0("draw.", dt$draw)]

# let's now calculate the statistics across all 100 draws
dt <- rbind(dt, dt[, .(pfpr = mean(pfpr), draw = "draw.mean"), by = .(areaId)])
dt <- rbind(dt, dt[, .(pfpr = sd(pfpr), draw = "draw.sd"), by = .(areaId)])

dt <- dcast(dt, areaId ~ draw, value.var = "pfpr" )
#View(dt[areaId %in% unique(areasf$id)])

# ##Make some maps ####
# area.data = merge(areasf, dt, by.x = "id", by.y = "areaId", all=FALSE)
# plot.data <- area.data[order(area.data$order), ]
# p1 = ggplot(data = plot.data, aes(x=long, y=lat, group = group))
# p2 = p1 + geom_polygon(data = bioko, aes(x = long, y = lat, group = group), color = "black", fill="grey", size = 0.25)
# map <- p2 + geom_polygon(data = plot.data, aes(x = long, y = lat, group = group, fill = draw.sd), color = NA, size = 0.25) +
#   scale_fill_gradient(name="pfpr", low="yellow", high="red", limits=c(0,.15)) +
#   geom_polygon(data = bioko, aes(x = long, y = lat, group = group), color = "black", fill=NA, size = 0.25) +
#   theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(), panel.background=element_blank(), legend.position=c(0.2, 0.8))
# map


# Store data ####
fwrite(dt, file = paste0(here("data/clean/"), "pfpr_draws.csv"))

# Code Graveyard ####
# and for good measure, let's compare to the mean + sd surfaces ...

# filename <- paste0(here("data/raw/Output_popMCD/"), "Prevalence_mean.tif")
# dat <- raster::raster(filename)
# dat.1 <- focal(dat, w = matrix(1, nc = 3, nr = 3), fun = mean, na.rm = TRUE)
# 
# pfpr.dt = data.table(
#   #coordinates(centroids), # lat-long coordinates of the centroids
#   areaId = centroids$FID, # this is the same as areaId
#   pfpr = extract(dat, centroids), # extract out the values of the raster map at the centroids of each pixel
#   pfpr.edges.1 = extract(dat.1, centroids) # extract out values of the edge-corrected raster map at centroids
# )
# # unfortunately, we are missing some data around the edges when we apply the extraction
# # and try to find the value of the raster at all of the centroids;
# # instead, we build out the edges using the focal function() and sample from there
# # not elegant, but here's how we fill in the missing edges on the same surface as before
# pfpr.dt[is.na(pfpr) & !is.na(pfpr.edges.1)]$pfpr <- pfpr.dt[is.na(pfpr) & !is.na(pfpr.edges.1)]$pfpr.edges.1
# pfpr.dt$pfpr.edges.1 <- NULL
# 
# dt$mn <- pfpr.dt$pfpr
# 
# filename <- paste0(here("data/raw/Output_popMCD/"), "Prevalence_sd.tif")
# dat <- raster::raster(filename)
# dat.1 <- focal(dat, w = matrix(1, nc = 3, nr = 3), fun = mean, na.rm = TRUE)
# 
# pfpr.dt = data.table(
#   #coordinates(centroids), # lat-long coordinates of the centroids
#   areaId = centroids$FID, # this is the same as areaId
#   pfpr = extract(dat, centroids), # extract out the values of the raster map at the centroids of each pixel
#   pfpr.edges.1 = extract(dat.1, centroids) # extract out values of the edge-corrected raster map at centroids
# )
# # unfortunately, we are missing some data around the edges when we apply the extraction
# # and try to find the value of the raster at all of the centroids;
# # instead, we build out the edges using the focal function() and sample from there
# # not elegant, but here's how we fill in the missing edges on the same surface as before
# pfpr.dt[is.na(pfpr) & !is.na(pfpr.edges.1)]$pfpr <- pfpr.dt[is.na(pfpr) & !is.na(pfpr.edges.1)]$pfpr.edges.1
# pfpr.dt$pfpr.edges.1 <- NULL
# 
# dt$sd <- pfpr.dt$pfpr
# 
# # Looks like we're pretty close...?
# dt[, .(areaId, draw.mean, draw.sd, mn, sd)][areaId %in% unique(areasf$id)]