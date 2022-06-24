rm(list=ls());gc()

# Load Libraries ####
library(data.table)
library(Matrix)
library(MASS)
library(here)

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

# care seeking behavior
FeverPf = 0.1116336
TreatPf = 0.602

# PfPR data
pfpr.data <- fread(here("data/pfpr_draws.csv"))
pfpr.data <- merge(pfpr.data, pop.data[year == 2018, .(map.area)], by = "map.area", all = FALSE)

# Dynamical Model Parameters
f = 0.3
q = 0.9
a = f*q
b = 0.55
c = 0.15
r = 1/200 # rate at which people become cured
eta = 1/32 # rate at which prophylaxis wears off
p = 0.9 # fraction of surviving mosquitoes
g = 1 - p # fraction of dying mosquitoes
eip <- 11
peip = p^eip # fraction of mosquitoes who survive incubation period
rho = FeverPf*TreatPf # probability of clearing infection through treatment cascade

# Define Time Spent at Risk Matrix
# Used for calibrating PfPR across the model
#
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

# Derive EIR ####
# Set random seed for drawing the surface from the joint posterior distribution #
# pfpr.draw.ix <- sample(1:100, 1)
# # Draw from PfPR surface
# pfpr.draw.colname = paste0("draw.", pfpr.draw.ix)
# pfpr.data.draw = pfpr.data[,c("map.area", pfpr.draw.colname), with = FALSE]
# colnames(pfpr.data.draw)[2] = "pfpr"
# pop.data <- merge(pop.data[year == 2018], pfpr.data.draw, by = "map.area")
# pfpr.input <- c(pop.data$pfpr, 0.43)

# use mean pfpr
pfpr.input <- c(pfpr.data$draw.mean, 0.43)

# See Ruktanonchai et al. (2016) and Supplementary Information for derivation
odds.vector <- r/(1-rho)*pfpr.input/(1-(1+rho*r/eta/(1-rho))*pfpr.input)
# Force of Infection (FOI) or "happenings rate" h
h.FOI <- MASS::ginv(TaR.matrix) %*% odds.vector
h.FOI[which(h.FOI < 0)] <- 0

# Calculate Mosquito Parameters ####
# we begin by calculating the fraction of infected mosquitoes, the sporozoite rate, from the Ross-Macdonald equations + Kappa
# Total population, including visitors
H.visitors <- t(c(pop.data[year == 2018]$pop, 0) %*% TaR.matrix)
# Sick population, including visitors
X.visitors <- t((c(pop.data[year == 2018]$pop, 0) * pfpr.input)  %*% TaR.matrix)
# This is the number of people who are sick, including visitors and residents both
kappa <- X.visitors/H.visitors
kappa <- as.vector(kappa)
# Define kappa off-island, based on PfPR off-island
kappa[242] <- .43
# z.spz <- peip*a*c*kappa/(p*a*c*kappa + (1-p))
z.spz <- peip*a*c*kappa/(a*c*kappa + (1-p))
# this is Z/M, but we currently do not know M
M = h.FOI*H.visitors/a/b/z.spz
M[242] = 0 # for off-island

Z = z.spz * M
Z[242] = 0 # for off-island
Y = Z/peip
Y[242] = 0 # for off-island

M <- as.vector(M)
Y <- as.vector(Y)
Z <- as.vector(Z)

# Derive Lambda, set the emergence rate of mosquitoes in each of the patches
# this is Lambda, calculated based on equilibrium value for the emergence process
# Lambda = M*(1-p)/p
lambda <- M * (1-p)
lambda[242] = 0 # for off-island
lambda <- as.vector(lambda)



n.patch <- 242 # 241 + 1 : the last patch is off-island

# set the EIR off-island
eg.eir <- h.FOI[242]/b


# Mosquito Parameters ####
# After setting up the patches, we populate them with mosquitoes 
# and define the emergence rate of mosquitoes in each patch
psi <- diag(1, nrow = n.patch)

# mosy_pars <- mosquito_rm_conpars(N = n.patch,
#                                  lambda = lambda.matrix,
#                                  psi = psi,
#                                  EIP = rep(11,365),
#                                  M = M,
#                                  Y = Y,
#                                  Z = Z)

# Human Parameters ####
# After setting up the patches, we populate them with people
patch.human.pop <- c(pop.data[year == 2018]$pop, 0) 

# total number of humans
n.humans <- sum(patch.human.pop)


# messing around




# # mosquito parameters
# f <- 0.3
# q <- 0.9
# eip <- 14
# EIP <- eip + 1
# 
# p <- 0.9
# g <- 1 - p
# peip <- p^EIP
# 
# # human parameters
# b <- 0.55
# c <- 0.15
# 
# # transmission parameters
# kappac <- kappa * c
# 
# h <- as.vector(h.FOI)
# N <- as.vector(H.visitors)
# 
# # equilibrium solutions
# # Z_new <- (h*N) / a # from FOI equation h = fqZ/N
# Z_new <- rep(50, 242)
# Z_new[242] <- 0
# Y_new <- Z_new / peip # exact for continuous or discrete time equations
# Y_new[242] <- 0
# M_new <- (Z_new*(g + (f*q*kappac))) / (f*q*kappac*peip)
# M_new[242] <- 0
# lambda_new <- g*M_new
# lambda_new[242] = 0

# --------------------------------------------------------------------------------
#   MicroMoB: stochastic
# --------------------------------------------------------------------------------

tmax <- 365 * 3

mod <- make_MicroMoB(tmax = tmax, p = n.patch)
setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = eip - 1, p = p, psi = psi, M = M, Y = Y, Z = Z)
setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE)

# setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = eip, p = p, psi = psi, M = M_new, Y = Y_new, Z = Z_new)
# setup_aqua_trace(model = mod, lambda = lambda_new, stochastic = FALSE)

simout <- data.table::CJ(day = 1:tmax, state = c('M', 'Y', 'Z'), patch = 1:n.patch, value = NaN)
simout <- simout[c('M', 'Y', 'Z'), on="state"]
data.table::setkey(simout, day)

# run it
pb <- progress_bar$new(total = tmax)

while(get_tnow(mod) <= tmax) {
  
  mod$mosquito$kappa <- kappa * c
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

# abs(M - simout[day == tmax & state == 'M', value])
