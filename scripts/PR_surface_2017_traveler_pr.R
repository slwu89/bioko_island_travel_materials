library(data.table)
library(sp)
library(raster)
library(rgdal)
library(here)
library(ggplot2)
here()

# load in centroids
centroids <- rgdal::readOGR(here("data/raw/mapArea_centroids/mapareagrid_centroids.shp"))
# load in prevalence conditioning on recent travel
dat <- raster::raster("data/raw/Outputs_prevalence_traveler_2017/t_Prevalence_mean.tif")
# do not attempt to reconcile the two by setting the crs to be the same
# somehow, because the extents are so different, this really does mess up a lot
# but if I don't do anything before calling the extract

# Use the focal() function to deal with the edge effects
# https://gis.stackexchange.com/questions/187410/how-do-i-get-rid-of-edge-effects-while-using-focal-in-r-to-smooth-a-raster
dat.1 <- focal(dat, w = matrix(1, nc = 3, nr = 3), fun = mean, na.rm = TRUE)

pfpr.dt = data.table(
  coordinates(centroids), # lat-long coordinates of the centroids
  areaId = centroids$OBJECTID, # this is the same as areaId
  pfpr = extract(dat, centroids), # extract out the values of the raster map at the centroids of each pixel
  pfpr.edges.1 = extract(dat.1, centroids) # extract out values of the edge-corrected raster map at centroids
)
# unfortunately, we are missing some data around the edges when we apply the extraction
# and try to find the value of the raster at all of the centroids;
# instead, we build out the edges using the focal function() and sample from there
# not elegant, but here's how we fill in the missing edges on the same surface as before
pfpr.dt[is.na(pfpr) & !is.na(pfpr.edges.1)]$pfpr <- pfpr.dt[is.na(pfpr) & !is.na(pfpr.edges.1)]$pfpr.edges.1
pfpr.dt$pfpr.edges.1 <- NULL

fwrite(pfpr.dt, here("data/clean/Mean_Traveler_PR_2017"))


# What the plots end up looking like:
areas.inh <- merge(bi.centroids, pop.summary, by.x = "OBJECTID", by.y = "areaId", all = FALSE)
areas.inh <- merge(areas.inh, pfpr.dt[,.(areaId, pfpr.travel = pfpr)], by.x = "OBJECTID", by.y = "areaId")

ggplot() + 
  geom_polygon(data = bi.outline, aes(x = long, y = lat, group = group), color = "black", fill="grey", size = 0.25) + 
  geom_point(data = data.frame(areas.inh), mapping = aes(x = coords.x1, y = coords.x2, color = 100*pfpr.travel), 
             shape = 15, size = 2.1) + 
  scale_color_gradient(name="PFPR | Travel(%)",
                       low="yellow", high="red", limits=c(0, 50)) +
  # move the scale bar
  coord_fixed() + # aspect ratio
  theme(aspect.ratio = 1, 
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(size =18),
        legend.text = element_text(size = 14),
        legend.position = c(.25, .8))
