####
# Clean Data - Regression Model of Trip Frequency
#
# Load in the aggregated trip data from clean/aggregated_2015_2018_travel_data.csv
#
# Model trip frequency based on probability of leaving: trip.counts/n
# and combine with the duration of the study period 56 days
# Return: a vector, broken down by map.area, of how frequently people leave home
#
# September 3, 2019
# Last Updated September 27, 2020
#
####

library(data.table)
library(here)
library(MASS)
library(nnet)

# Load aggregated travel data ####
travel.data <- fread(here("data/aggregated_2015_2018_travel_data.csv"))
# merge with travel distance data set:
reg.travel.dist <- fread(here("data/travel_dist_by_region.csv"))
travel.data <- merge(travel.data, reg.travel.dist, by = "map.area")
# Format:
# map.area
# ad2
# n - number surveyed, by year
# trip.counts - number of people who report leaving (for 2018, this is different from the sum over trips)
# t_eg etc. - number of reported trips to each destination
# pop - population, 2015 and 2018 censuses
# year - year of survey
# 

# here's syntax for how to aggregate over years
# travel.data[,.(n = sum(n), trip.counts = sum(trip.counts), t_eg = sum(t_eg), ti_mal = sum(ti_mal)), by = .(map.area, ad2)]

# list of all unique map.areas
map.area.list <- sort(unique(travel.data$map.area))

# Part 1: Modeling travel frequency for each map.area ####

# Frequency and Probability of Leaving ####
# Model: trip.counts/n
# The probability of leaving, modeled as a binomial choice
# Use as covariates the population, admin2 region, and time to Malabo:
# there are a few places where the total number of trips is more than the number of people
# another way to model this would be with a Poisson-type model of the frequency of leaving
# for now, let's just find the probability of leaving at all
freq.model.data <- travel.data
freq.model.data[trip.counts > n]$trip.counts <- freq.model.data[trip.counts > n]$n

h <- glm( cbind(trip.counts, n - trip.counts) ~ pop + ad2 + dist.mal,
          data = freq.model.data,
          family = binomial(link = logit))
freq.model.data$leave.prob <- h$fitted.values
freq.model.data$leave.freq <- h$fitted.values/56

# how different are the results from if we had aggregated over the years?
# how different are the results from if we had used trip.counts = sum for the 2018 data?
# makes almost no difference

freq.model.data.clean <- freq.model.data[, c("map.area", "year", "leave.prob"), with = FALSE]

fwrite(freq.model.data.clean, file = here("data/clean/trip_frequency_model_estimates.csv"))


# Perform draws
# sampled.model <- mvrnorm(100, coefficients(h), vcov(h))
# freq.samples.dt <- matrix(0, nrow = length(freq.model.data$map.area), ncol = 103)
# freq.samples.dt[,1] <- freq.model.data$map.area # list of map.areas, across all 4 years
# freq.samples.dt[,2] <- freq.model.data$year
# freq.samples.dt[,3] <- h$fitted.values # mean value
# h.draw <- h # here's where we alter the coefficients and resample
# for(i in 1:100){
#   h.draw$coefficients <- sampled.model[i,]
#   freq.samples.dt[,(i+3)] <- predict(h.draw, freq.model.data, type = "response")
# }
# freq.samples.dt <- as.data.table(freq.samples.dt)
# colnames(freq.samples.dt) <- c("map.area", "year", "draw.mean", paste0("draw.",as.character(1:100)))

# Write out frequency model estimates ####
# fwrite(freq.samples.dt, file = here("data/clean/trip_frequency_model_estimates.csv"))
