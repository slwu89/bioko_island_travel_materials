####
# Bioko Island EIR Surface calculator
# Requires aggregating a lot of scripts that I've already created
# Updated from the version used for ASTMH, this time accommodating the 
# disambiguated TaR matrix.
#
# This is going to be practice for taking in and compiling all of the data that 
# we already have, using it to form a single TaR matrix and a single lambda surface
#
# August 8, 2019 - in preparation for writing the spatial uncertainty manuscript
#
####

# Load Libraries ####
library(here)
library(data.table, quietly = TRUE) # The workhorse of this exercise
library(pscl, quietly = TRUE) # Used for fitting the travel model
library(boot, quietly = TRUE)
library(MASS, quietly = TRUE)
library(nnet, quietly = TRUE)
library(rgdal, quietly = TRUE)

# Read in and format the travel data ####
# Read in population denominator
denom <- data.table(read.csv("/Users/dtcitron/Documents/MASH/MASH-Main/MASH-dev/DanielCitron/Bioko_Island_Cluster_Simulations/BI_sim_setup_data/bioko_areas.csv"))
# PfPR Estimates - from Carlos's paper using Su's MAP model
BI <- as.data.table(read.csv("/Users/dtcitron/Documents/MASH/Bioko_Macro/Travel_Model_Fitting/Importations - Carlos Paper/eta.csv"))
# Travel Survey Data
travel <- data.table(read.csv("/Users/dtcitron/Documents/MASH/MASH-Main/MASH-dev/DanielCitron/Bioko_Island_Cluster_Simulations/BI_sim_setup_data/summaries.csv"))
# Travel Times - proxy for distance between things
times <- data.table(read.csv("/Users/dtcitron/Documents/MASH/MASH-Main/MASH-dev/DanielCitron/Bioko_Island_Cluster_Simulations/BI_sim_setup_data/travel_times.csv"))
# rename some of the columns, and reorder
times <- times[, 5:199]
colnames(times)[2:195] <- sapply(colnames(times)[2:195], FUN = function(s){
  return(substring(s, 2, 5))}
)
setcolorder(times, c(1,order(as.integer(colnames(times)[2:195]))+1))
times <- times[order(areaId)]


# Visualizations set-up ####
# When we want to make maps, need to read in and set up area-level grid
#bioko<-readOGR("/Users/dtcitron/Documents/MASH/MASH-Main/MASH-dev/DanielCitron/Bioko_Island_Cluster_Simulations/BI_maps/bioko", "bioko_admin1")
#areas_inh<-readOGR("/Users/dtcitron/Documents/MASH/MASH-Main/MASH-dev/DanielCitron/Bioko_Island_Cluster_Simulations/BI_maps/areas_inh", "areas_inh")
#areasf<-fortify(areas_inh, region = "areaId")


# Merge and combine relevant data ####
# Combine travel with pfprsu based on areaId
travel.model.data <- merge(travel, BI[, c("areaId", "to_map", "prt_map", "prall_map", "te_map", "pre_map")], by = "areaId")
# old version
#travel <- merge(travel, pfpr[,.(areaId, pfpr)], by = "areaId")

# Combine geographic data with travel based on areaId
travel.model.data<- merge(travel.model.data, denom[,c("areaId", "pop", "lon", "lat")], by = "areaId")

# Define a vector of AreaIDs ####
# We'll be using this to find the indexes corresponding to the 
areaIds <- sort(unique(travel.model.data$areaId))

# Subsetting data from the raw data loaded in previously
travel.model.data <- travel.model.data[, .(areaId, lon, lat, X, Y, # areaID and location
                                           # A2 region and whether or not the patch is near Malabo
                                           ad2, malabo,
                                           # Population, PfPR
                                           pop, prall_map, pre_map, prt_map,
                                           # Travel Frequency
                                           te_map, to_map,
                                           # Cluster ID
                                           clusId,
                                           # Number surveyed, number surveyed with Pf+
                                           n, pf,
                                           # Number surveyed who reported fever, got treatment
                                           pffv, pftto,
                                           # Number of times traveled to each destination, 2015-2017
                                           to, ti_ban, ti_mal, ti_lub, ti_ria, ti_mok, ti_ure
)]
# Replace a LOT of NAs with other values:
travel.model.data[is.na(travel.model.data$n),]$n <- 0
travel.model.data[is.na(to),]$to <- 0
travel.model.data[is.na(ti_ban),]$ti_ban <- 0
travel.model.data[is.na(ti_lub),]$ti_lub <- 0
travel.model.data[is.na(ti_mal),]$ti_mal <- 0
travel.model.data[is.na(ti_mok),]$ti_mok <- 0
travel.model.data[is.na(ti_ria),]$ti_ria <- 0
travel.model.data[is.na(ti_ure),]$ti_ure <- 0

# Relabel all peri-urban (Punto Europa) areas in Malabo as "Peri"
travel.model.data[ad2 == "Malabo" & malabo == "non-Malabo"]$ad2 <- "Peri"
travel.model.data$ad2 <- as.factor(travel.model.data$ad2)