####
# Raw data - travel times data from 2015-2017
# Compute mean travel times to each ad2 location
#
# NB: we are going to use distances for now, because I don't think it makes much of a difference?
#
# Output data table
#
# September 3, 2019
#
####

library(here)
library(data.table)
library(rgdal)


# Read in data ####
# centroid locations, used for calculating distances
centroids <- rgdal::readOGR(here("data/raw/mapArea_centroids/mapareagrid_centroids.shp"))
centroids <- as.data.table(centroids)
# travel data information; includes areaId and ad2 for 2015-2018
travel.data <- fread(here("data/clean/aggregated_2015_2018_travel_data.csv"))
# merge them together
# needed to fix this, because otherwise we leave out 2514 ...
centroids <- merge(travel.data[year == 2018 | (year == 2015 & areaId == 2514),.(areaId, ad2)], centroids, by.x = "areaId", by.y = "OBJECTID", all = FALSE)


# Fill out matrix of travel distances (not times)
areaId.list <- centroids$areaId
trip.distances <- matrix(0, nrow = length(areaId.list), ncol = (length(areaId.list)+1))
trip.distances[, 1] <- areaId.list
# loop over rows
for(i in 1:length(areaId.list)){
  # loop over columns
  for(j in 1:(length(areaId.list))){
    x1 = centroids[areaId == areaId.list[[i]]]$coords.x1
    y1 = centroids[areaId == areaId.list[[i]]]$coords.x2
    x2 = centroids[areaId == areaId.list[[j]]]$coords.x1
    y2 = centroids[areaId == areaId.list[[j]]]$coords.x2
    d = sqrt((x1 - x2)^2 + (y1 - y2)^2)
    trip.distances[i,j+1] <- d
  }
}

colnames(trip.distances) <- c("areaId", as.character(areaId.list))

# goal is to write down travel distance to each region:
# find all areaIds in each target region:
ban.areaIds <- travel.data[ad2 == "Baney"]$areaId
ban.ixs <- sapply(ban.areaIds, function(y) which(y == areaId.list))
mal.areaIds <- travel.data[ad2 %in% c("Peri","Malabo")]$areaId
mal.ixs <- sapply(mal.areaIds, function(y) which(y == areaId.list))
lub.areaIds <- travel.data[ad2 == "Luba"]$areaId
lub.ixs <- sapply(lub.areaIds, function(y) which(y == areaId.list))
ria.areaIds <- travel.data[ad2 == "Riaba"]$areaId
ria.ixs <- unlist(sapply(ria.areaIds, function(y) which(y == areaId.list)))
ure.areaIds <- travel.data[ad2 == "Ureka"]$areaId
ure.ixs <- sapply(ure.areaIds, function(y) which(y == areaId.list))
mok.areaIds <- travel.data[ad2 == "Moka"]$areaId
mok.ixs <- sapply(mok.areaIds, function(y) which(y == areaId.list))

# calculate average distance between 152 and each of the 6 locations ####
reg.distances <- matrix(0, ncol = 7, nrow = length(areaId.list))
for(i in 1:length(areaId.list)){
  reg.distances[i,] <- c(areaId.list[i], 
                         mean(trip.distances[i, ban.ixs+1]),
                         mean(trip.distances[i, lub.ixs+1]),
                         mean(trip.distances[i, mal.ixs+1]),
                         mean(trip.distances[i, mok.ixs+1]),
                         mean(trip.distances[i, ria.ixs+1]),
                         mean(trip.distances[i, ure.ixs+1]))
}
reg.distances <- as.data.table(reg.distances)
colnames(reg.distances) <- c("areaId", "dist.ban", "dist.lub", "dist.mal", "dist.mok", "dist.ria", "dist.ure")

# add off-island travel distances ####
# distance to malabo + geographical distance to mainland eg
reg.distances$dist.eg <- reg.distances$dist.mal + 240000

# Merge and Save ####
fwrite(reg.distances, file = here("data/clean/travel_dist_by_region.csv"))
trip.distances.dt <- as.data.table(trip.distances)
fwrite(trip.distances.dt, file = here("data/clean/travel_dist_by_areaId.csv"))


