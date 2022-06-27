rm(list=ls());gc()

# Load Libraries ####
library(data.table)
library(Matrix)
library(MASS)
library(here)
library(MicroMoB)
library(progress)
library(ggplot2)

# Load Data and Set Parameter Values ####
# population data
pop.data <- fread(here("data/aggregated_2015_2018_travel_data.csv"))
# travel frequencies
trip.freq.data <- fread(here("data/trip_frequency_model_estimates.csv"))
# travel destination selection model
trip.dest.data <- fread(here("data/negative_binomial_predictions_by_destination_region.csv")) 

# travel duration
trip.duration.eg <- 1/0.04713604 # rate of return from mainland eg to bioko
trip.duration.bi <- 1/0.09670302 # rate of return from trips on bioko

# PfPR data
pfpr.data <- fread(here("data/pfpr_draws.csv"))
pfpr.data <- merge(pfpr.data, pop.data[year == 2018, .(map.area)], by = "map.area", all = FALSE)

# get the region - to - pixel matrix
source(here("scripts/region_to_maparea_mapping.R"))
# convert the trip.dest.data into a matrix
trip.dest.mat <- as.matrix(trip.dest.data[year == 2018, .(t_eg, ti_ban, ti_lub, ti_mal, ti_mok, ti_ria, ti_ure)])
trip.dest.mat <- rbind(trip.dest.mat, matrix(c(1,0,0,0,0,0,0), ncol = 7))
# take the matrix product of trip.dest.mat with the region-to-pixel matrix toget pixel-to-pixel:
movement.matrix <- trip.dest.mat %*% reg.2.pixel

# Trip durations, by destination across the island and off-island
trip.duration <- c(rep(trip.duration.bi, 241), trip.duration.eg)

# Frequencies at which people leave home
# Divide by 56 days to transform the probability of leaving into 
# the frequency of leaving during the 8-week study period
trip.freq <- trip.freq.data[year == 2018, c("map.area", "leave.prob"), with = FALSE]
colnames(trip.freq)[2] <- "freq"
trip.freq <- c(trip.freq[order(map.area)]$freq/56, 1)

# Build full Time Spent at Risk matrix
TaR.matrix <- diag(1, nrow = 242, ncol = 242)
for (i in 1:242){
  TaR.matrix[i,] <- (movement.matrix[i,]*trip.duration)/(sum(movement.matrix[i,]*trip.duration) + 1/trip.freq[i])
  TaR.matrix[i,i] <- 1 - sum(TaR.matrix[i,])
  TaR.matrix[i,] <- TaR.matrix[i,]/sum(TaR.matrix[i,])
}
TaR.matrix <- TaR.matrix[1:241, 1:241]

pfpr <- pfpr.data$draw.mean

# mosquito parameters
f = 0.3
q = 0.9
a = f*q
p = 0.9 # mosquito survival
g = 1 - p # mosquito mortality
eip <- 11
peip = p^eip # fraction of mosquitoes who survive incubation period

# transmission parameters
b = 0.55
c = 0.15

# human infection parameters
FeverPf = 0.1116336
TreatPf = 0.602
r = 1/200 # rate at which people become cured
eta = 1/32 # rate at which prophylaxis wears off
rho = FeverPf*TreatPf # probability of clearing infection through treatment cascade

# See Ruktanonchai et al. (2016) and Supplementary Information for derivation
odds.vector <- r/(1-rho)*pfpr/(1-(1+rho*r/eta/(1-rho))*pfpr)

# Force of Infection (FOI) or "happenings rate" h
h.FOI <- MASS::ginv(TaR.matrix) %*% odds.vector
h.FOI[which(h.FOI < 0)] <- 0
h.FOI <- as.vector(h.FOI)

# Calculate Mosquito Parameters ####
# we begin by calculating the fraction of infected mosquitoes, the sporozoite rate, from the Ross-Macdonald equations + Kappa
# Total population, including visitors
H.visitors <- t(pop.data[year == 2018, pop] %*% TaR.matrix)
# Sick population, including visitors
X.visitors <- t((pop.data[year == 2018, pop] * pfpr)  %*% TaR.matrix)

# kappa: net infectiousness of human pop at each place
kappa <- X.visitors/H.visitors * c
kappa <- as.vector(kappa)

H.visitors <- as.vector(H.visitors)
n_patch <- length(kappa)


# --------------------------------------------------------------------------------
# mosy input parameters
# --------------------------------------------------------------------------------

psi <- diag(n_patch)
surv <- peip

Z <- (h.FOI * H.visitors) / (a * b)
Y <- Z / surv
M <- (Z*(g + (f*q*p*kappa))) / (f*q*kappa*p*surv) # M from kappa
lambda <- g*M

# run sim
tmax <- 365 * 2

mod <- make_MicroMoB(tmax = tmax, p = n_patch)
setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = eip, p = p, psi = psi, M = M, Y = Y, Z = Z)
setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE)

simout <- data.table::CJ(day = 1:tmax, state = c('M', 'Y', 'Z'), value = NaN, patch = 1:n_patch)
simout <- simout[c('M', 'Y', 'Z'), on="state"]
data.table::setkey(simout, day, patch)

# run it
pb <- progress_bar$new(total = tmax)

while(get_tnow(mod) <= tmax) {
  
  mod$mosquito$kappa <- kappa
  step_aqua(model = mod)
  step_mosquitoes(model = mod)
  
  simout[day == get_tnow(mod) & state == 'M', value := mod$mosquito$M]
  simout[day == get_tnow(mod) & state == 'Y', value := mod$mosquito$Y]
  simout[day == get_tnow(mod) & state == 'Z', value := mod$mosquito$Z]
  
  mod$global$tnow <- mod$global$tnow + 1L
  pb$tick()
}

simout_summary <- simout[, .(value = sum(value)), by = .(state, day)]

ggplot(simout_summary) +
  geom_line(aes(x = day, y = value, color = state)) +
  facet_wrap(. ~ state, scales = "free")


mod$mosquito$Z / surv
mod$mosquito$Y
