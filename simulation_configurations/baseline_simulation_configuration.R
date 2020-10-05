##
#
# Script for Outputting Baseline Bioko Island Simulations
#
# Daniel T Citron
#
##

rm(list=ls());gc()

# Load Libraries ####
library(data.table)
library(Matrix)
library(MASS)
library(here)

# Load macro.pfsi, library for simulations
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
# If not installed:
#devtools::install_github(repo = "https://github.com/dd-harp/MASH",subdir = "macro.pfsi")
library(macro.pfsi)

# Load Data and Set Parameter Values ####
# population data
pop.data <- fread(here("data/clean/aggregated_2015_2018_travel_data.csv"))
# travel frequencies
trip.freq.data <- fread(here("data/clean/trip_frequency_model_estimates.csv"))
# travel destination selection model
trip.dest.data <- fread(here("data/clean/negative_binomial_predictions_by_destination_region.csv")) 

# travel duration
trip.duration.eg <- 1/0.04713604 # rate of return from mainland eg to bioko
trip.duration.bi <- 1/0.09670302 # rate of return from trips on bioko

# care seeking behavior
FeverPf = 0.1116336
TreatPf = 0.602

# PfPR data
pfpr.data <- fread(here("data/clean/pfpr_draws.csv"))
pfpr.data <- merge(pfpr.data, pop.data[year == 2018, .(areaId)], by = "areaId", all = FALSE)


# Dynamical Model Parameters
a = 0.3*0.9
b = 0.55
c = 0.15
r = 1./200 # rate at which people become cured
eta = 1./32 # rate at which prophylaxis wears off
p = 0.9 # fraction of surviving mosquitoes
g = 1 - p # fraction of dying mosquitoes
peip = p^11 # fraction of mosquitoes who survive incubation period
rho = FeverPf*TreatPf # probability of clearing infection through treatment cascade


# Define Time Spent at Risk Matrix
# Used for calibrating PfPR across the model
#
# get the region - to - pixel matrix
source(here("scripts/region_to_areaId_mapping.R"))
# convert the trip.dest.data into a matrix
tar.draw.ix = paste0("draw.","mean")
trip.dest.mat <- as.matrix(trip.dest.data[year == 2018 & draw == tar.draw.ix, .(t_eg, ti_ban, ti_lub, ti_mal, ti_mok, ti_ria, ti_ure)])
trip.dest.mat <- rbind(trip.dest.mat, matrix(c(1,0,0,0,0,0,0), ncol = 7))
# take the matrix product of trip.dest.mat with the region-to-pixel matrix toget pixel-to-pixel:
movement.matrix <- trip.dest.mat %*% reg.2.pixel

# Trip durations, by destination across the island and off-island
trip.duration <- c(rep(trip.duration.bi, 241), trip.duration.eg)

# Frequencies at which people leave home
# Divide by 56 days to transform the probability of leaving into 
# the frequency of leaving during the 8-week study period
freq.draw.ix = paste0("draw.","mean")
trip.freq <- trip.freq.data[year == 2018, c("areaId", freq.draw.ix), with = FALSE]
colnames(trip.freq)[2] <- "freq"
trip.freq <- c(trip.freq[order(areaId)]$freq/56, 1)

# Build full Time Spent at Risk matrix
TaR.matrix <- diag(1, nrow = 242, ncol = 242)
for (i in 1:242){
  TaR.matrix[i,] <- (movement.matrix[i,]*trip.duration)/(sum(movement.matrix[i,]*trip.duration) + 1/trip.freq[i])
  TaR.matrix[i,i] <- 1 - sum(TaR.matrix[i,])
  TaR.matrix[i,] <- TaR.matrix[i,]/sum(TaR.matrix[i,])
}

# Derive EIR ####
# Set random seed for drawing the surface from the joint posterior distribution #
set.seed(1)
pfpr.draw.ix <- sample(1:100, 1)
# Draw from PfPR surface
pfpr.draw.colname = paste0("draw.", pfpr.draw.ix)
pfpr.data.draw = pfpr.data[,c("areaId", pfpr.draw.colname), with = FALSE]
colnames(pfpr.data.draw)[2] = "pfpr"
pop.data <- merge(pop.data[year == 2018], pfpr.data.draw, by = "areaId")
pfpr.input <- c(pop.data$pfpr, 0.43)

# See Ruktanonchai et al. (2016) and Supplementary Information for derivation
odds.vector <- r/(1-rho)*pfpr.input/(1-(1+rho*r/eta/(1-rho))*pfpr.input)
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
# Define kappa off-island, based on PfPR off-island
kappa[242] <- .43
z.spz <- peip*a*c*kappa/(p*a*c*kappa + (1-p))
# this is Z/M, but we currently do not know M
M = h.FOI*H.visitors/a/b/z.spz

Z = z.spz * M
Z[242] = 0 # for off-island
Y = Z/peip
Y[242] = 0 # for off-island

# Derive Lambda, set the emergence rate of mosquitoes in each of the patches
# this is Lambda, calculated based on equilibrium value for the emergence process
Lambda = M*(1-p)/p
Lambda[242] = 0 # for off-island

# PfSI Parameters ####
# Parameters definie probabilities of symptomatic malaria and treatment seeking
pfsi_pars <- pfsi_parameters(FeverPf = 0.1116336, TreatPf = 0.602)

# Patch Parameters ####
# set up patches (n is how many patches we have)
n.patch <- 242 # 241 + 1 : the last patch is off-island

# set the EIR off-island
eg.eir <- h.FOI[242]/b
patch_pars <- patches_parameters(move = movement.matrix,
                                 bWeightZoo = rep(0,n.patch),
                                 bWeightZootox = rep(0,n.patch),
                                 reservoir = c(rep(F,(n.patch-1)), T),
                                 res_EIR = c(rep(0,(n.patch-1)),eg.eir)
)

# Mosquito Parameters ####
# After setting up the patches, we populate them with mosquitoes 
# and define the emergence rate of mosquitoes in each patch
psi <- Matrix::sparseMatrix(i = {},j = {},x = 0.0,dims = c(n.patch,n.patch))
diag(psi) <- rep(1,n.patch)
lambda.matrix = t(matrix(Lambda, nrow = n.patch, ncol = 365))
mosy_pars <- mosquito_rm_conpars(N = n.patch,
                                 lambda = lambda.matrix,
                                 psi = psi,
                                 EIP = rep(11,365),
                                 M = M,
                                 Y = Y,
                                 Z = Z)


# Human Parameters ####
# After setting up the patches, we populate them with people
patch.human.pop <- c(pop.data[year == 2018]$pop, 0) 

# total number of humans
n.humans <- sum(patch.human.pop) 

# Randomly set initial infection status of each human host malaria prevalence in each patch
pfpr <- pfpr.input
set.seed(1)
init_state <- unlist(
  mapply(FUN = function(n,pr){
    sample(x = c("I","S"),size = n,replace = T,prob = c(pr,1-pr))
  },
  n=patch.human.pop,
  pr=pfpr.input,
  SIMPLIFY = F)
)

# Define Patch IDs, for where people go
# These IDs correspond to indices in a vector - need it to be 0-indexed for c++
patch_id <- rep(0:(n.patch-1), times=patch.human.pop)

# Set uniform biting weights
bweights <- rep(1, n.humans)

# Set mean trip durations based on destination
trip.durations <- c(rep(trip.duration.bi, 241), trip.duration.eg)
# Set trip frequencies - this is set according to one's home origin patch
trip.freqs <- trip.freqs <- rep(trip.freq[1:241], times = patch.human.pop[1:241])


# The data structure for constructing the human pop
human_pars <- vector("list", n.humans)
for(i in 1:n.humans){
  human_pars[[i]] <- human_pfsi_conpars(id = i-1,
                                        home_patch_id = patch_id[i],
                                        trip_duration = trip.durations,
                                        trip_frequency = trip.freqs[i],
                                        bweight = bweights[i],
                                        age = 20,
                                        state = init_state[i],
                                        bite_algorithm = 0)
}

# Run Simulation ####

# Define Output Path Names #
log_pars <- list()
h_inf <- here("data/simulation_outputs", paste0("pfsi_", 1, ".csv"))
log_pars[[1]] <- list(outfile = h_inf, 
                      key = "pfsi",
                      header = paste0(c("time",
                                        "patch",
                                        unlist(lapply(c("S","I","P"),function(x){
                                          paste0(x,c("_visitor","_resident_home","_resident_away"))}
                                        )),
                                        "incidence_resident",
                                        "incidence_traveller"),
                                      collapse = ",")
)

mosy <-  here("data/simulation_outputs", paste0("mosy_", 1, ".csv"))
log_pars[[2]] <- list(outfile = mosy,
                      key = "mosquito",
                      header = paste0(c("time",
                                        "state",
                                        paste0("patch",1:n.patch)),
                                      collapse = ","))

# Set random seed
set.seed(1)
# 
run_macro(tmax = 7*365,
          human_pars = human_pars,
          mosquito_pars = mosy_pars,
          patch_pars = patch_pars,
          model_pars = pfsi_pars,
          log_streams = log_pars,
          vaxx_events = NULL,
          verbose = T)


# Analyzing the output ####
library(ggplot2)

h_inf <- here("data/simulation_outputs", paste0("pfsi_", 1, ".csv"))
dt <- fread(h_inf)

# Create a new data table to merge onto the simulation output, to track population denominators over time
areaId.list <- sort(pop.data[year == 2018]$areaId)
pop.dt <- data.table(patch = c(0:(241-1)), areaId = areaId.list)
pop.dt <- merge(pop.dt, pop.data[year == 2018, .(areaId, pop)], by = "areaId")
# merge, to use the pop column as a denominator when calculating fractions susceptible, infected, protected
dt <- merge(dt, pop.dt, by = "patch")
dt[, s := (S_resident_home + S_resident_away)/pop, by = c("time" , "patch" , "time")]
dt[, i := (I_resident_home + I_resident_away)/pop, by = c("time" , "patch" , "time")]
dt[, p := (P_resident_home + P_resident_away)/pop, by = c("time" , "patch" , "time")]

h <- melt(dt[areaId %in%  c(152, 207, 220,335,502, 644, 1175,2199,2457)],
          id.vars = c("time", "areaId"),
          measure.vars = c("s","i","p"),
          value.name = "fraction")

# Create a plot
ggplot(data = h) +
  geom_point(mapping = aes(x = time, y = fraction, color = variable), shape = 20, size = .01) +
  facet_wrap(~areaId)
