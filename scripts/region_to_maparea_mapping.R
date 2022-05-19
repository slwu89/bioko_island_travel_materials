####
# Map Admin 2 "destination regions" from the travel survey data
# onto map.areas within those destination regions
#
# We need a matrix to perform this mapping; this script produces that matrix.
# The matrix will be weighted by map.area populations (from 2018), but we can
# convert the matrix to an unweighted matrix.
#
# The matrix should have 241 + 1 columns, where each of the 241 columns 
# corresponds to an map.area (from 2018) and the last column corresponds to off-island travel
# The matrix should have 7 rows, each of the 7 rows corresponds to a destination region,
# including off-island travel.
# The value of each entry for the ith map.area and the jth destination region should be
# the population-weighted probability of going to map.area i given that one has gone to
# destination region j.
#
# September 13, 2019
# Updated October 29, 2019
#
####

library(data.table)
library(here)
# load data ####
travel.data <- fread(here("data/aggregated_2015_2018_travel_data.csv"))
# subset data on 2018
dat <- travel.data[year == 2018,.(map.area, ad2, pop)]


# list of all map.areas in 2018 ####
map.area.list <- dat$map.area
# determine fractions of population in each location, where the denominator is the regional population ####
dat$pop.frac <- 0
dat[ad2 == "Baney"]$pop.frac <- dat[ad2 == "Baney"]$pop/sum(dat[ad2 == "Baney"]$pop)
dat[ad2 == "Luba"]$pop.frac <- dat[ad2 == "Luba"]$pop/sum(dat[ad2 == "Luba"]$pop)
dat[ad2 %in% c("Peri","Malabo")]$pop.frac <- dat[ad2 %in% c("Peri","Malabo")]$pop/sum(dat[ad2 %in% c("Peri","Malabo")]$pop)
dat[ad2 == "Moka"]$pop.frac <- dat[ad2 == "Moka"]$pop/sum(dat[ad2 == "Moka"]$pop)
dat[ad2 == "Riaba"]$pop.frac <- dat[ad2 == "Riaba"]$pop/sum(dat[ad2 == "Riaba"]$pop)
dat[ad2 == "Ureka"]$pop.frac <- dat[ad2 == "Ureka"]$pop/sum(dat[ad2 == "Ureka"]$pop)
# 
# Construct the matrix ####
reg.2.pixel <- matrix(0, nrow = 7, ncol = (length(map.area.list)+1))
# Fill in the matrix
# the order of indexes for destinations is: Off-Island, Baney, Luba, Malabo, Moka, Riaba, Ureka
off.ixs <- c(242) # this is the map.area index for off-island
reg.2.pixel[1,off.ixs] <- 1
ban.ixs <- match(dat[ad2 == "Baney"]$map.area, map.area.list)
reg.2.pixel[2,ban.ixs] <- dat[ad2 == "Baney"]$pop.frac
lub.ixs <- match(dat[ad2 == "Luba"]$map.area, map.area.list)
reg.2.pixel[3,lub.ixs] <- dat[ad2 == "Luba"]$pop.frac
mal.ixs <- match(dat[ad2 %in% c("Peri","Malabo")]$map.area, map.area.list)
reg.2.pixel[4,mal.ixs] <- dat[ad2 %in% c("Peri","Malabo")]$pop.frac
mok.ixs <- match(dat[ad2 == "Moka"]$map.area, map.area.list)
reg.2.pixel[5,mok.ixs] <- dat[ad2 == "Moka"]$pop.frac
ria.ixs <- match(dat[ad2 == "Riaba"]$map.area, map.area.list)
reg.2.pixel[6,ria.ixs] <- dat[ad2 == "Riaba"]$pop.frac
ure.ixs <- match(dat[ad2 == "Ureka"]$map.area, map.area.list)
reg.2.pixel[7,ure.ixs] <- dat[ad2 == "Ureka"]$pop.frac
# 
# # Check:
# rowSums(reg.2.pixel)
# sum(colSums(reg.2.pixel) == 0)
