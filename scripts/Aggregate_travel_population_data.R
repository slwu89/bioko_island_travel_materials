####
# Raw data - Aggregate population and travel data by pixel + year
#
# Load in the summary of the 2018 MIS travel data
# Load in the 2015-2017 travel data as well
# Combine the data, but keep each year separated!
# 
# The problem to solve is that there are some differences between the 2018 and 2015 census
# We can account for the fact that some places are inhabited in one data set and
# uninhabited in another data set by
# 
# We will assume that the census data is the same for 2015-2017, and then for 2018
# Alternatively, we might assume that population trends linearly from year to year
# (which is simple but probably isn't true).
#
# September 9, 2019
#
####


# Load libraries ####
library(data.table)
library(here)

# Data table formatting ####
# n is the number of people sampled
# trip.counts is the number of people who report leaving home <- we use this to calculate travel frequency
# t_eg is off-island travel; ti_ban is on-island travel to Baney, etc <- we use this to calculate destination choice
# pop is population
# year is the year of the data

# Load 2015-2017 data ####
travel.data.2015.2017 <- fread(here("data/raw/2015-2017_survey_data/summaries.csv"))
# Adjust ad2 categorization
pixels.peri <- c(152, 153, 207, 209, 211, 212, 218, 219, 220, 270, 271,
                 329, 330, 386, 387, 443, 445, 447, 448, 502, 503, 504, 
                 505, 506, 507, 560, 564, 571, 573 ,574, 617, 618, 630, 
                 633, 634, 676, 677, 693, 734, 735, 736, 792, 793, 794, 
                 851, 910, 969, 970, 1027, 1028)
travel.data.2015.2017[areaId %in% pixels.peri]$ad2 <- "Peri"
travel.data.2015.2017[areaId == 397]$ad2 <- "Baney"
travel.data.2015.2017[areaId %in% c(218,219,220)]$ad2 <- "Malabo"
# Population
pop.by.area.2015.2017 <- fread(here("data/raw/2015-2017_survey_data/bioko_areas.csv"))
pop.by.area.2015.2017 <- pop.by.area.2015.2017[,.(areaId, pop)]


# 2015 data subset ####
travel.data.2015 <- travel.data.2015.2017[,c("areaId", "ad2", "n15", "to15", "ti_ban_15", "ti_lub_15", "ti_mal_15", "ti_mok_15", "ti_ria_15", "ti_ure_15")]
colnames(travel.data.2015) <- c("areaId", "ad2", "n", "t_eg", "ti_ban", "ti_lub", "ti_mal", "ti_mok", "ti_ria", "ti_ure")
travel.data.2015[is.na(n)]$ti_ure <- 0
travel.data.2015[is.na(n)]$ti_ria <- 0
travel.data.2015[is.na(n)]$ti_mok <- 0
travel.data.2015[is.na(n)]$ti_mal <- 0
travel.data.2015[is.na(n)]$ti_lub <- 0
travel.data.2015[is.na(n)]$ti_ban <- 0
travel.data.2015[is.na(n)]$t_eg <- 0
travel.data.2015[is.na(n)]$n <- 0
travel.data.2015[, trip.counts := sum(t_eg, ti_ban, ti_lub, ti_mal, ti_mok, ti_ria, ti_ure), by = "areaId"]
travel.data.2015 <- merge(travel.data.2015, pop.by.area.2015.2017, by = "areaId")
travel.data.2015$year <- 2015
setcolorder(travel.data.2015, c("areaId", "ad2", "n", "trip.counts", "t_eg", "ti_ban", "ti_lub", "ti_mal", "ti_mok", "ti_ria", "ti_ure", "pop", "year"))
# 2016 data subset ####
travel.data.2016 <- travel.data.2015.2017[,c("areaId", "ad2", "n16", "to16", "ti_ban_16", "ti_lub_16", "ti_mal_16", "ti_mok_16", "ti_ria_16", "ti_ure_16")]
colnames(travel.data.2016) <- c("areaId", "ad2", "n", "t_eg", "ti_ban", "ti_lub", "ti_mal", "ti_mok", "ti_ria", "ti_ure")
travel.data.2016[is.na(n)]$ti_ure <- 0
travel.data.2016[is.na(n)]$ti_ria <- 0
travel.data.2016[is.na(n)]$ti_mok <- 0
travel.data.2016[is.na(n)]$ti_mal <- 0
travel.data.2016[is.na(n)]$ti_lub <- 0
travel.data.2016[is.na(n)]$ti_ban <- 0
travel.data.2016[is.na(n)]$t_eg <- 0
travel.data.2016[is.na(n)]$n <- 0
travel.data.2016[, trip.counts := sum(t_eg, ti_ban, ti_lub, ti_mal, ti_mok, ti_ria, ti_ure), by = "areaId"]
travel.data.2016 <- merge(travel.data.2016, pop.by.area.2015.2017, by = "areaId")
travel.data.2016$year <- 2016
setcolorder(travel.data.2016, c("areaId", "ad2", "n", "trip.counts", "t_eg", "ti_ban", "ti_lub", "ti_mal", "ti_mok", "ti_ria", "ti_ure", "pop", "year"))
# 2017 data subset ####
travel.data.2017 <- travel.data.2015.2017[,c("areaId", "ad2", "n17", "to17", "ti_ban_17", "ti_lub_17", "ti_mal_17", "ti_mok_17", "ti_ria_17", "ti_ure_17")]
colnames(travel.data.2017) <- c("areaId", "ad2", "n", "t_eg", "ti_ban", "ti_lub", "ti_mal", "ti_mok", "ti_ria", "ti_ure")
travel.data.2017[is.na(n)]$ti_ure <- 0
travel.data.2017[is.na(n)]$ti_ria <- 0
travel.data.2017[is.na(n)]$ti_mok <- 0
travel.data.2017[is.na(n)]$ti_mal <- 0
travel.data.2017[is.na(n)]$ti_lub <- 0
travel.data.2017[is.na(n)]$ti_ban <- 0
travel.data.2017[is.na(n)]$t_eg <- 0
travel.data.2017[is.na(n)]$n <- 0
travel.data.2017[, trip.counts := sum(t_eg, ti_ban, ti_lub, ti_mal, ti_mok, ti_ria, ti_ure), by = "areaId"]
travel.data.2017 <- merge(travel.data.2017, pop.by.area.2015.2017, by = "areaId")
travel.data.2017$year <- 2017
setcolorder(travel.data.2017, c("areaId", "ad2", "n", "trip.counts", "t_eg", "ti_ban", "ti_lub", "ti_mal", "ti_mok", "ti_ria", "ti_ure", "pop", "year"))

# 2018 data ####
mis.raw <- fread(here("data/raw/2018_travel_data/mis2018travel.csv"))
# Open data set for catching trip counts data ####
travel.data.2018 <- data.table(areaId = sort(unique(mis.raw$areaId)))
# Count the total number of people interviewed in each areaId in 2018:
travel.data.2018 <- merge(travel.data.2018, mis.raw[, .(n = .N), by = areaId], by = "areaId", all = TRUE)

# For each areaId, find the number of people who left home
travel.data.2018 <- merge(travel.data.2018, mis.raw[travelledisland == "yes" | travelledEG == "yes", .(trip.counts = .N), by = areaId, ],
                          by = "areaId", all= TRUE)
travel.data.2018[is.na(trip.counts)]$trip.counts <- 0
# 
# # For each areaID, find the number of people who left home to go elsewhere on the island
# travel.data.2018 <- merge(travel.data.2018, mis.raw[travelledisland == "yes", .(counts.bi = .N), by = areaId, ],
#                           by = "areaId", all= TRUE)
# travel.data.2018[is.na(counts.bi)]$counts.bi <- 0
# 
# For each areaID, find the number of people who left home to go to the mainland
travel.data.2018 <- merge(travel.data.2018, mis.raw[travelledEG == "yes", .(counts.eg = .N), by = areaId, ],
                          by = "areaId", all= TRUE)
travel.data.2018[is.na(counts.eg)]$counts.eg <- 0


# Read in the off-island travel data ####
# The problem here is that we have missing information for off-island data
# so for the most part there are fewer detailed trips to the mainland
# than there are people who report going to the mainland!
# Discounting multiple mainland trips, we'll just use the reported
# off-island travel question from mis.raw
# eg.trips <- fread(here("data/raw/2018_travel_data/trips_to_EG.csv"))
# holder <- eg.trips[, .(off.island.trips = .N), by = "areaId"][order(areaId)]
# travel.data.2018 <- merge(travel.data.2018, holder, by = "areaId", all = TRUE)

# Read in the on-island travel data ####
bi.trips <- fread(here("data/raw/2018_travel_data/trips_on_BI.csv"))
# now we need to map the ad4 units onto ad2 units
# first read in ad4 shapefile
library(sp)
library(raster)
library(maptools)
library(rgdal)
ad4 <- rgdal::readOGR(here("data/raw/admin4shp/admin4V19.shp"))
ad4.2.ad2 <- as.data.table(ad4[,c("admin2", "admin2ID", "admin4ID")])
colnames(ad4.2.ad2) <- c("admin2", "admin2Id", "admin4Id")
# now merge to include ad2 destinations
bi.trips <- merge(bi.trips, ad4.2.ad2, by = "admin4Id")
# create a separate column for destination region
bi.trips$dest.reg <- bi.trips$admin2
# now we designate certain destination ad4's as Moka and Ureka
# Ureka's ad4 ids: L260 ("Ureka") - 2 surveyed
# mis.raw[community == "Ureka", c("community", "admin4ID")
# Moka's ad4 ids: L184 ("Moka Malabo") - 77 surveyed
# mis.raw[community == "Moka Malabo", c("community", "admin4ID")]
# L266 (Moka Bioko) - 231 surveyed
# mis.raw[community == "Moka Bioko", c("community", "admin4ID")]
bi.trips[admin4Id == "L260"]$dest.reg <- "Ureka"
bi.trips[admin4Id %in% c("L184", "L266")]$dest.reg <- "Moka"

# Transform data: count number of trips for each areaId to each destination region
holder <- bi.trips[, c("areaId", "dest.reg")]
holder <- dcast(holder, areaId ~ dest.reg, length)

# Put everything together
travel.data.2018 <- merge(travel.data.2018, holder, by = "areaId", all = TRUE)
travel.data.2018[is.na(travel.data.2018)] <- 0

# Add 2018 populations
pop.by.area <- fread(here("data/clean/pop_areas2018.csv"))
# need to be careful to include all inhabited pixels
travel.data.2018 <-  merge(travel.data.2018, pop.by.area[!is.na(pop) & pop > 0,.(areaId, pop)], by = "areaId", all= TRUE)
travel.data.2018 <- travel.data.2018[, lapply(.SD, function(x){ifelse(is.na(x), 0, x)})]
travel.data.2018[areaId == 704]$pop <- 5 # for some reason the population here is zeroed out in the census
travel.data.2018$year <- 2018

# Note that for 2018 there is a weird discrepancy between trip.counts and the sum over 
# the individual reported trips to all destinations.  This is because the trip.counts
# counts people who left home to travel on-island + off-island as a single trip
# We are effectively measuring the probability of leaving home > 0 times.
# We can alternatively replace trip.counts with the sum total of detailed trip counts.

travel.data.2018 <- merge(travel.data.2018, travel.data.2015.2017[,.(areaId, ad2)], by = "areaId", all.x = TRUE, all.y = FALSE)
travel.data.2018[areaId %in% c(341, 400, 410, 456, 469)]$ad2 <- "Baney"
pixels.peri.2018 <- c(327, 328, 385, 388, 389, 446, 567, 568, 570, 572, 519, 619, 626, 627, 687, 795, 909, 968, 1261)
travel.data.2018[areaId %in% pixels.peri.2018]$ad2 <- "Peri"
pixels.luba <- c(1732, 2027, 2084, 2085, 2131, 2135, 2200, 2204, 2307, 2377, 2497)
travel.data.2018[areaId %in% pixels.luba]$ad2 <- "Luba"
pixels.riaba <- c(2043, 2110, 2227, 2286, 2335, 2337, 2338, 2341, 2393, 2394, 2397, 2403, 2404, 2753)
travel.data.2018[areaId %in% pixels.riaba]$ad2 <- "Riaba"
travel.data.2018[areaId == 3500]$ad2 <- "Ureka"

# rename and reorder columns:
colnames(travel.data.2018) <- c("areaId", "n", "trip.counts", "t_eg", "ti_ban", "ti_lub", "ti_mal", "ti_ria", "ti_ure", "ti_mok", "pop", "year", "ad2")
setcolorder(travel.data.2018, c("areaId", "ad2", "n", "trip.counts", "t_eg", "ti_ban", "ti_lub", "ti_mal", "ti_mok", "ti_ria", "ti_ure", "pop", "year"))

# Plotting, showing where each pixel falls on the map 
# areas_inh <- raster::shapefile(here("data/raw/mapArea_centroids/mapareagrid_centroids.shp")) #
# areas_inh <- areas_inh[areas_inh$OBJECTID %in% travel.data.2018$areaId,]
# h <- data.frame(areas_inh)
# h2 <- travel.data.2018[, .(OBJECTID = areaId, ad2)]
# h <- merge(h, h2, by = "OBJECTID" , all = TRUE)
# #h[h$OBJECTID == 2497,]$ad2 <- "hippo"
# ggplot() +
#   geom_point(data = h, mapping = aes(x = coords.x1, y = coords.x2, color = ad2))

# Concatenate these data sets:
aggregated.travel.data <- rbind(travel.data.2015, travel.data.2016, travel.data.2017, travel.data.2018)
fwrite(aggregated.travel.data, here("data/clean/aggregated_2015_2018_travel_data.csv"))
