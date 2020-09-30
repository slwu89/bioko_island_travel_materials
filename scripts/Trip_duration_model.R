####
# Clean Data - Regression Model of Trip Duration
#
# Load in the trip duration data from raw/2018_travel_data/trips_on_BI.csv
# and raw/2018_travel_data/trips_to_EG.csv
#
# Model trip duration as a function of destination for each ad2 destination
#
# August 12, 2019
#
####

library(here)
library(data.table)

# Fit data for trips to EG ####

# Load the data:
dat <- fread(here("data/clean/BI_travel_duration_data.csv"))
# split data into on-island and off-island
eg.dat <- dat[destination == "on_island"]
bi.dat <- dat[destination == "off_island" & !is.na(nights)]

# Fit to an exponential decay, using maximum likelihood:
f <- optimise(f = function(l){
  sum(dexp(x = eg.dat[!is.na(nights)]$nights, rate = l, log = TRUE))
  }, 
  interval = c(0,1), 
  maximum = TRUE)
# This is the decay rate that maximizes the likelihood
lambda.eg <- f$maximum # 0.04713604

# Make a histogram
h <- hist(eg.dat$nights, breaks = seq(0,380,7))
# and the fit itself
x <- seq(0,370,1)
y <- dexp(x = x, rate = 0.0471359, log = FALSE)
# Plot the data, just to show that it works okay:
plot(h$breaks[1:54] + 3.5, h$density, col = "blue")
points(x,y, add =TRUE, axes = FALSE)


# QWerage duration for travel anywhere on the island:

f.bi.all <- optimise(f = function(l){
    sum(dexp(x = bi.dat$nights, rate = l, log = TRUE))
  }, 
  interval = c(0,1), 
  maximum = TRUE)
lambda.bi.all <- f.bi.all$maximum # 0.09670302
1/f.bi.all$maximum # average is a 10-day stay

# make a histogram, and plot everything again
h <- hist(bi.dat$nights, breaks = seq(0,271,1))
x <- seq(0,270,1)
y <- dexp(x = x, rate = lambda.bi.all, log = FALSE)
# Plot the data, just to show that it works okay:
plot(h$breaks[1:271] + 1, h$density, col = "blue")
points(x,y, add =TRUE, axes = FALSE)
