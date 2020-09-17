##
#
# Script for Outputting Low-importations Bioko Island Simulations
# For ASTMH
# Propagating through spatial uncertainty 
#
# Daniel T Citron
# 10/31/19
#
##
#
# A single run of the simulation
# Designed to be called using a qsub job on the cluster
# Pass in a random seed, which will serve as the seed for the simulation.
# The random seed will also serve as a seed for drawing randomly from 1 to 100 surfaces.
#
##

# This is for Low Importations, what happens when we decrease exposure to Malaria risk on the mainland
# We are going to use the real initial conditions, for making beautiful plots!
#
# To clarify: half as much malaria risk on the mainland
#
# We still calibrate M Z Y Lambda based on the old FOI, calculated with high off-island exposure
# But we will change the off-island FOI when parameterizing off-island EIR;
# ALso when initializing PfPR IC's

rm(list=ls());gc()

seed <- as.integer(commandArgs()[8])
print(seed)
# Load Libraries ####
library(data.table)
library(ggplot2)
library(Matrix)
library(MASS)

setwd("/ihme/malaria_modeling/dtcitron")
library(here, lib.loc = "/ihme/malaria_modeling/dtcitron/Rlibs")

library(macro.pfsi, lib.loc = here("Rlibs"))

# Create dummy output, for testing the qsubber script:####
# dt <- data.table(x = c(1,2,3,4), y = c(2,4,6,8))
# fwrite(dt, here("SpatialUncertainty/workflow_test/qsub_outputs",paste0("qsub_test_",seed,".csv")))

# Set random seed for drawing the surface from the joint posterior distribution ####
set.seed(seed)
pfpr.draw.ix <- sample(1:100, 1)


# Load Data ####
# population data
pop.data <- fread(here("SpatialUncertainty/data_clean/aggregated_2015_2018_travel_data.csv"))
# travel frequency
trip.freq.data <- fread(here("SpatialUncertainty/data_clean/trip_frequency_model_estimates.csv"))
# trip.freq.data[year == 2018, .(areaId, year, draw.mean)]

# travel destination selection model
trip.dest.data <- fread(here("SpatialUncertainty/data_clean/negative_binomial_predictions_by_destination_region.csv")) 
# trip.dest.data[year == 2018 & draw == "draw.mean"]

# travel duration
trip.duration.eg <- 1/0.04713604 # rate of return from mainland eg to bioko
trip.duration.bi <- 1/0.09670302 # rate of return from trips on bioko

# care seeking behavior
FeverPf = 0.1116336
TreatPf = 0.602

# PfPR data
pfpr.data <- fread(here("SpatialUncertainty/data_clean/pfpr_draws.csv"))
pfpr.data <- merge(pfpr.data, pop.data[year == 2018, .(areaId)], by = "areaId", all = FALSE)


# Dynamical Model Parameters ####
a = 0.3*0.9
b = 0.55
c = 0.15
r = 1./200 # rate at which people become cured
eta = 1./32 # rate at which prophylaxis wears off
p = 0.9 # fraction of surviving mosquitoes
g = 1 - p # fraction of dying mosquitoes
peip = p^11 # fraction of mosquitoes who survive incubation period
rho = FeverPf*TreatPf # probability of clearing infection through treatment cascade


# TaR Matrix ####
# get the region - to - pixel matrix
source(here("SpatialUncertainty/data_clean/region_to_areaId_mapping.R"))
# convert the trip.dest.data into a matrix
tar.draw.ix = paste0("draw.","mean")
trip.dest.mat <- as.matrix(trip.dest.data[year == 2018 & draw == tar.draw.ix, .(t_eg, ti_ban, ti_lub, ti_mal, ti_mok, ti_ria, ti_ure)])
trip.dest.mat <- rbind(trip.dest.mat, matrix(c(1,0,0,0,0,0,0), ncol = 7))
# take the matrix product of trip.dest.mat with the region-to-pixel matrix toget pixel-to-pixel:
movement.matrix <- trip.dest.mat %*% reg.2.pixel

# NB: the movement matrix is the probability of going from each areaId to each other areaId;
# we need to combine it with the other travel data to obtain the full TaR matrix
TaR.matrix <- diag(1, nrow = 242, ncol = 242)
# vector of trip durations across the island and off-island
trip.duration <- c(rep(trip.duration.bi, 241), trip.duration.eg)
# vector of frequencies at which people leave home
# where we divide by 56 to transform the probability of leaving into 
# the frequency of leaving during the study period

freq.draw.ix = paste0("draw.","mean")
trip.freq <- trip.freq.data[year == 2018, c("areaId", freq.draw.ix), with = FALSE]
colnames(trip.freq)[2] <- "freq"
trip.freq <- c(trip.freq[order(areaId)]$freq/56, 1)
for (i in 1:242){
  TaR.matrix[i,] <- (movement.matrix[i,]*trip.duration)/(sum(movement.matrix[i,]*trip.duration) + 1/trip.freq[i])
  TaR.matrix[i,i] <- 1 - sum(TaR.matrix[i,])
  TaR.matrix[i,] <- TaR.matrix[i,]/sum(TaR.matrix[i,])
}


# Derive EIR ####
pfpr.draw.colname = paste0("draw.", pfpr.draw.ix)
pfpr.data.draw = pfpr.data[,c("areaId", pfpr.draw.colname), with = FALSE]
colnames(pfpr.data.draw)[2] = "pfpr"
pop.data <- merge(pop.data[year == 2018], pfpr.data.draw, by = "areaId")
pfpr.input <- c(pop.data$pfpr, .43) # Low Importations: This is the PfPR input; it is not the same as kappa, which we do use to parameterize lambda etc

odds.vector <- r/(1-rho)*pfpr.input/(1-(1+rho*r/eta/(1-rho))*pfpr.input)
h.FOI <- MASS::ginv(TaR.matrix) %*% odds.vector
h.FOI[which(h.FOI < 0)] <- 0


# Derive Initial Conditions for Mosquitoes ####
# we begin by calculating the fraction of infected mosquitoes, the sporozoite rate, from the Ross-Macdonald equations + Kappa
# Total population, including visitors
H.visitors <- t(c(pop.data[year == 2018]$pop, 0) %*% TaR.matrix)
# Sick population, including visitors
X.visitors <- t((c(pop.data[year == 2018]$pop, 0) * pfpr.input)  %*% TaR.matrix)
# This is the number of people who are sick, including visitors and residents both
kappa <- X.visitors/H.visitors
kappa[242] <- .43/2 # Low importations - this is where we turn down the off-island kappa
z.spz <- peip*a*c*kappa/(p*a*c*kappa + (1-p))
# this is Z/M, but we currently do not know M
M = c(h.FOI[1:241]*H.visitors[1:241]/a/b/z.spz[1:241], 0)  # Local Residual

Z = z.spz * M
Z[242] = 0 # for off-island
Y = Z/peip
Y[242] = 0 # for off-island


# Derive Lambda ####
# We use this quantity to set the emergence rate of mosquitoes in each of the patches
# this is Lambda, calculated based on equilibrium value for the emergence process
Lambda = M*(1-p)/p
Lambda[242] = 0 # for off-island


# PfSI Parameters ####
# vector of parameters 
pfsi_pars <- pfsi_parameters(FeverPf = 0.1116336, TreatPf = 0.602)


# Patch Parameters ####
# set up patches (n is how many patches we have)
n.patch <- 242 # 241 + 1 : the last patch is off-island

# set the EIR off-island
eg.eir <- r/(1-rho)*(.43/2)/(1-(1+rho*r/eta/(1-rho))*(.43/2))/b # Low importations - this should be about half of the actual value #h.FOI[242]/b
patch_pars <- patches_parameters(move = movement.matrix,
                                 bWeightZoo = rep(0,n.patch),
                                 bWeightZootox = rep(0,n.patch),
                                 reservoir = c(rep(F,(n.patch-1)), T),
                                 res_EIR = c(rep(0,(n.patch-1)),eg.eir)
)


# Mosy Parameters ####
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
# number of people in each patch
patch.human.pop <- c(pop.data[year == 2018]$pop, 0) 
# malaria prevalence in each patch
pfpr <- pfpr.input # true ICs
# total number of humans
n.humans <- sum(patch.human.pop) 


# sample S or I for each person
# set seed first, for the sake of reproducibility - this is the same as the seed passed to this script
set.seed(seed)
init_state <- unlist(
  mapply(FUN = function(n,pr){
    sample(x = c("I","S"),size = n,replace = T,prob = c(pr,1-pr))
  },
  n=patch.human.pop,
  pr=pfpr,
  SIMPLIFY = F)
)

# Define Patch IDs, for where people go
# These IDs correspond to indices in a vector - need it to be 0-indexed for c++
patch_id <- rep(0:(n.patch-1), times=patch.human.pop)

# Set uniform biting weights; this could follow any density on the positive reals (gamma, log-normal, weibull, etc.)
bweights <- rep(1,n.humans)

# Set mean trip durations based on destination
trip.durations <- c(rep(trip.duration.bi, 241), trip.duration.eg)
# Set trip frequencies - this is set according to one's home origin patch
trip.freqs <- trip.freqs <- rep(trip.freq[1:241], times = patch.human.pop[1:241])


# the data structure that will be passed down to C++ to construct the human pop
human_pars <- vector("list",n.humans)
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


# Define Path Names ####
log_pars <- list()
#h_inf <- paste0(path,"pfsi.csv")
h_inf <- here("SpatialUncertainty/ASTMH19/lowimport_trueICs/sim_output", paste0("pfsi_", seed, ".csv"))
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
#mosy <- paste0(path,"mosy.csv")
mosy <-  here("SpatialUncertainty/ASTMH19/lowimport_trueICs/sim_output", paste0("mosy_", seed, ".csv"))
log_pars[[2]] <- list(outfile = mosy,
                      key = "mosquito",
                      header = paste0(c("time",
                                        "state",
                                        paste0("patch",1:n.patch)),
                                      collapse = ","))


# Run Simulation ####
# Set random seed, for reproducibility
# Use the seed passed to this program
set.seed(seed)
# 
run_macro(tmax = 10*365,
          human_pars = human_pars,
          mosquito_pars = mosy_pars,
          patch_pars = patch_pars,
          model_pars = pfsi_pars,
          log_streams = log_pars,
          vaxx_events = NULL,
          verbose = T)


# Analyzing the output - comment out later ####

# h_inf <- here("SpatialUncertainty/ASTMH19/lowimport_trueICs/sim_output", paste0("pfsi_", 2, ".csv"))
# dt <- fread(h_inf)
# 
# # Create a new data table to merge onto the simulation output, to track population denominators over time
# areaId.list <- sort(pop.data[year == 2018]$areaId)
# pop.dt <- data.table(patch = c(0:(241-1)), areaId = areaId.list)
# pop.dt <- merge(pop.dt, pop.data[year == 2018, .(areaId, pop)], by = "areaId")
# # merge, to use the pop column as a denominator when calculating fractions susceptible, infected, protected
# dt <- merge(dt, pop.dt, by = "patch")
# dt[, s := (S_resident_home + S_resident_away)/pop, by = c("time" , "patch" , "time")]
# dt[, i := (I_resident_home + I_resident_away)/pop, by = c("time" , "patch" , "time")]
# dt[, p := (P_resident_home + P_resident_away)/pop, by = c("time" , "patch" , "time")]
# 
# h <- melt(dt[areaId %in%  c(335, 2694)],#152, 207, 220,335,502, 644, 1175,2199,2457)],
#           id.vars = c("time", "areaId"),
#           measure.vars = c("s","i","p"),
#           value.name = "fraction")
# 
# ggplot(data = h) +
#   geom_point(mapping = aes(x = time, y = fraction, color = variable), shape = 20, size = .01) +
#   facet_wrap(~areaId) +
#   ylim(0,.25)

