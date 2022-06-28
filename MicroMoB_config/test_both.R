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

pfpr <- c(pfpr.data$draw.mean, 0.43)

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
H.visitors <- t(c(pop.data[year == 2018]$pop, 0) %*% TaR.matrix)
# Sick population, including visitors
X.visitors <- t((c(pop.data[year == 2018]$pop, 0) * pfpr)  %*% TaR.matrix)

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


# --------------------------------------------------------------------------------
# human input parameters
# --------------------------------------------------------------------------------

patch.human.pop <- c(pop.data[year == 2018]$pop, 0)

# # total number of humans
n.humans <- sum(patch.human.pop)

# set the EIR off-island
eg.eir <- h.FOI[242]/b



# --------------------------------------------------------------------------------
#   FOI stuff
# --------------------------------------------------------------------------------

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


# initial conditions
human_state <- matrix(data = 0, nrow = n_patch, ncol = 3, dimnames = list(NULL, c('S', 'I', 'P')))
human_state[, 'I'] <- patch.human.pop * pfpr
human_state[, 'S'] <- patch.human.pop - human_state[, 'I']

tmax <- 365 * 2
mod <- make_MicroMoB(tmax = tmax, p = n_patch)

setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE)
setup_mosquito_RM(model = mod, stochastic = FALSE, f = f, q = q, eip = eip, p = p, psi = psi, nu = 25, M = M, Y = Y, Z = Z)
setup_humans_SIP(model = mod, stochastic = FALSE, theta = TaR.matrix, SIP = human_state, b = b, c = c, r = r, rho = rho, eta = eta)
setup_alternative_trace(model = mod)
setup_visitor_trace(model = mod)


n <- mod$global$n
p <- mod$global$p

# human quantities
H <- compute_H(mod)
x <- compute_x(mod)
wf <- compute_wf(mod)
Psi <- compute_Psi(mod)
W <- as.vector(t(Psi) %*% (wf * H))

# ok
sum(Psi == TaR.matrix)

# ok
W == H.visitors

# biting distribution matrix (n x p)
beta <- diag(wf, nrow = n, ncol = n) %*% Psi %*% diag(1/W, nrow = p, ncol = p)

# feeding rate
f <- compute_f(mod, B = B)

# human blood feeding fraction
q <- compute_q(mod, W = W, Wd = Wd, B = B)

# density of infective mosquitoes
Z <- compute_Z(mod)

# calculate EIR and kappa (mosy->human, human->mosy
EIR <- beta %*% (f*q*Z)
EIR <- as.vector(EIR)
h <- EIR * b


# --------------------------------------------------------------------------------
#   deterministic simulation
# --------------------------------------------------------------------------------

# initial conditions
human_state <- matrix(data = 0, nrow = n_patch, ncol = 3, dimnames = list(NULL, c('S', 'I', 'P')))
human_state[, 'I'] <- patch.human.pop * pfpr
human_state[, 'S'] <- patch.human.pop - human_state[, 'I']

tmax <- 365 * 2
mod <- make_MicroMoB(tmax = tmax, p = n_patch)

setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE)
setup_mosquito_RM(model = mod, stochastic = FALSE, f = f, q = q, eip = eip, p = p, psi = psi, nu = 25, M = M, Y = Y, Z = Z)
setup_humans_SIP(model = mod, stochastic = FALSE, theta = TaR.matrix, SIP = human_state, b = b, c = c, r = r, rho = rho, eta = eta)
setup_alternative_trace(model = mod)
setup_visitor_trace(model = mod)

# human output table
human_out <- data.table::CJ(day = 1:tmax, state = c('S', 'I', 'P'), patch = 1:n_patch, value = NaN)
human_out <- human_out[c('S', 'I', 'P'), on="state"]
data.table::setkey(human_out, day, patch)

# mosy output table
mosy_out <- data.table::CJ(day = 1:tmax, state = c('M', 'Y', 'Z'), patch = 1:n_patch, value = NaN)
mosy_out <- mosy_out[c('M', 'Y', 'Z'), on="state"]
data.table::setkey(mosy_out, day, patch)

pb <- progress_bar$new(total = tmax)

# run it
while (get_tnow(mod) <= tmax) {
  
  compute_bloodmeal(model = mod)
  
  mod$human$EIR[242] <- eg.eir
  
  step_aqua(model = mod)
  step_mosquitoes(model = mod)
  step_humans(model = mod)
  
  human_out[day==get_tnow(mod) & state=='S', value := mod$human$SIP[, 'S']]
  human_out[day==get_tnow(mod) & state=='I', value := mod$human$SIP[, 'I']]
  human_out[day==get_tnow(mod) & state=='P', value := mod$human$SIP[, 'P']]
  
  mosy_out[day == get_tnow(mod) & state == 'M', value := mod$mosquito$M]
  mosy_out[day == get_tnow(mod) & state == 'Y', value := mod$mosquito$Y]
  mosy_out[day == get_tnow(mod) & state == 'Z', value := mod$mosquito$Z]
  
  mod$global$tnow <- mod$global$tnow + 1L
  pb$tick()
}

human_summarized <- human_out[, .(value = sum(value)), by = .(state, day)]

ggplot(data = human_summarized) +
  geom_line(aes(x = day, y = value, color = state)) +
  facet_wrap(. ~ state, scales = "free") +
  ggtitle("deterministic humans")

mosy_out <- mosy_out[, .(value = sum(value)), by = .(state, day)]

ggplot(mosy_out) +
  geom_line(aes(x = day, y = value, color = state)) +
  facet_wrap(. ~ state, scales = "free") +
  ggtitle("deterministic mosy")









# rm(list=ls());gc()
# 
# # Load Libraries ####
# library(data.table)
# library(Matrix)
# library(MASS)
# library(here)
# 
# # Load Data and Set Parameter Values ####
# # population data
# pop.data <- fread(here("data/aggregated_2015_2018_travel_data.csv"))
# # travel frequencies
# trip.freq.data <- fread(here("data/trip_frequency_model_estimates.csv"))
# # travel destination selection model
# trip.dest.data <- fread(here("data/negative_binomial_predictions_by_destination_region.csv")) 
# 
# # travel duration
# trip.duration.eg <- 1/0.04713604 # rate of return from mainland eg to bioko
# trip.duration.bi <- 1/0.09670302 # rate of return from trips on bioko
# 
# # care seeking behavior
# FeverPf = 0.1116336
# TreatPf = 0.602
# 
# # PfPR data
# pfpr.data <- fread(here("data/pfpr_draws.csv"))
# pfpr.data <- merge(pfpr.data, pop.data[year == 2018, .(map.area)], by = "map.area", all = FALSE)
# 
# # Dynamical Model Parameters
# f = 0.3
# q = 0.9
# a = f*q
# b = 0.55
# c = 0.15
# r = 1/200 # rate at which people become cured
# eta = 1/32 # rate at which prophylaxis wears off
# p = 0.9 # fraction of surviving mosquitoes
# g = 1 - p # fraction of dying mosquitoes
# eip <- 11
# peip = p^eip # fraction of mosquitoes who survive incubation period
# rho = FeverPf*TreatPf # probability of clearing infection through treatment cascade
# 
# # Define Time Spent at Risk Matrix
# # Used for calibrating PfPR across the model
# #
# # get the region - to - pixel matrix
# source(here("scripts/region_to_maparea_mapping.R"))
# # convert the trip.dest.data into a matrix
# trip.dest.mat <- as.matrix(trip.dest.data[year == 2018, .(t_eg, ti_ban, ti_lub, ti_mal, ti_mok, ti_ria, ti_ure)])
# trip.dest.mat <- rbind(trip.dest.mat, matrix(c(1,0,0,0,0,0,0), ncol = 7))
# # take the matrix product of trip.dest.mat with the region-to-pixel matrix toget pixel-to-pixel:
# movement.matrix <- trip.dest.mat %*% reg.2.pixel
# 
# # Trip durations, by destination across the island and off-island
# trip.duration <- c(rep(trip.duration.bi, 241), trip.duration.eg)
# 
# # Frequencies at which people leave home
# # Divide by 56 days to transform the probability of leaving into 
# # the frequency of leaving during the 8-week study period
# trip.freq <- trip.freq.data[year == 2018, c("map.area", "leave.prob"), with = FALSE]
# colnames(trip.freq)[2] <- "freq"
# trip.freq <- c(trip.freq[order(map.area)]$freq/56, 1)
# 
# # Build full Time Spent at Risk matrix
# TaR.matrix <- diag(1, nrow = 242, ncol = 242)
# for (i in 1:242){
#   TaR.matrix[i,] <- (movement.matrix[i,]*trip.duration)/(sum(movement.matrix[i,]*trip.duration) + 1/trip.freq[i])
#   TaR.matrix[i,i] <- 1 - sum(TaR.matrix[i,])
#   TaR.matrix[i,] <- TaR.matrix[i,]/sum(TaR.matrix[i,])
# }
# 
# # Derive EIR ####
# # Set random seed for drawing the surface from the joint posterior distribution #
# # pfpr.draw.ix <- sample(1:100, 1)
# # # Draw from PfPR surface
# # pfpr.draw.colname = paste0("draw.", pfpr.draw.ix)
# # pfpr.data.draw = pfpr.data[,c("map.area", pfpr.draw.colname), with = FALSE]
# # colnames(pfpr.data.draw)[2] = "pfpr"
# # pop.data <- merge(pop.data[year == 2018], pfpr.data.draw, by = "map.area")
# # pfpr.input <- c(pop.data$pfpr, 0.43)
# 
# # use mean pfpr
# pfpr.input <- c(pfpr.data$draw.mean, 0.43)
# 
# # See Ruktanonchai et al. (2016) and Supplementary Information for derivation
# odds.vector <- r/(1-rho)*pfpr.input/(1-(1+rho*r/eta/(1-rho))*pfpr.input)
# # Force of Infection (FOI) or "happenings rate" h
# h.FOI <- MASS::ginv(TaR.matrix) %*% odds.vector
# h.FOI[which(h.FOI < 0)] <- 0
# 
# # Calculate Mosquito Parameters ####
# # we begin by calculating the fraction of infected mosquitoes, the sporozoite rate, from the Ross-Macdonald equations + Kappa
# # Total population, including visitors
# H.visitors <- t(c(pop.data[year == 2018]$pop, 0) %*% TaR.matrix)
# # Sick population, including visitors
# X.visitors <- t((c(pop.data[year == 2018]$pop, 0) * pfpr.input)  %*% TaR.matrix)
# # This is the number of people who are sick, including visitors and residents both
# kappa <- X.visitors/H.visitors
# # Define kappa off-island, based on PfPR off-island
# kappa[242] <- .43
# z.spz <- peip*a*c*kappa/(p*a*c*kappa + (1-p))
# # this is Z/M, but we currently do not know M
# M = h.FOI*H.visitors/a/b/z.spz
# M[242] = 0 # for off-island
# 
# Z = z.spz * M
# Z[242] = 0 # for off-island
# Y = Z/peip
# Y[242] = 0 # for off-island
# 
# M <- as.vector(M)
# Y <- as.vector(Y)
# Z <- as.vector(Z)
# 
# # Derive Lambda, set the emergence rate of mosquitoes in each of the patches
# # this is Lambda, calculated based on equilibrium value for the emergence process
# Lambda = M*(1-p)/p
# Lambda[242] = 0 # for off-island
# Lambda <- as.vector(Lambda)
# 
# 
# 
# n.patch <- 242 # 241 + 1 : the last patch is off-island
# 
# # set the EIR off-island
# eg.eir <- h.FOI[242]/b
# 
# 
# # Mosquito Parameters ####
# # After setting up the patches, we populate them with mosquitoes 
# # and define the emergence rate of mosquitoes in each patch
# psi <- diag(1, nrow = n.patch)
# 
# # mosy_pars <- mosquito_rm_conpars(N = n.patch,
# #                                  lambda = lambda.matrix,
# #                                  psi = psi,
# #                                  EIP = rep(11,365),
# #                                  M = M,
# #                                  Y = Y,
# #                                  Z = Z)
# 
# 
# # Human Parameters ####
# # After setting up the patches, we populate them with people
# patch.human.pop <- c(pop.data[year == 2018]$pop, 0) 
# 
# # total number of humans
# n.humans <- sum(patch.human.pop)