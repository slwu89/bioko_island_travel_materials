rm(list=ls());gc()
library(MicroMoB)
library(progress)
library(ggplot2)
library(deSolve)
library(data.table)
library(here)


# --------------------------------------------------------------------------------
# test SIP equilibrium for non spatial model
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


# --------------------------------------------------------------------------------
# testing FOI in space
# --------------------------------------------------------------------------------

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
TaR <- diag(1, nrow = 242, ncol = 242)
for (i in 1:242){
  TaR[i,] <- (movement.matrix[i,]*trip.duration)/(sum(movement.matrix[i,]*trip.duration) + 1/trip.freq[i])
  TaR[i,i] <- 1 - sum(TaR[i,])
  TaR[i,] <- TaR[i,]/sum(TaR[i,])
}

pfpr <- c(pfpr.data$draw.mean, 0.43)
N <- c(pop.data[year == 2018]$pop, 0)
# hack to fix problematically low pfpr values
pfpr[c(71,72,109,215,216)] <- pfpr[c(71,72,109,215,216)] * 1.5

# start of eq calc
odds <- r/(1-rho)*pfpr/(1-(1+rho*r/eta/(1-rho))*pfpr)

# Force of Infection (FOI) in each patch (not on each strata)
h_patch <- as.vector(MASS::ginv(TaR) %*% odds)

# EIR in patches
H <- as.vector(N %*% TaR)
Z <- (h_patch * H) / (a * b)

I <- pfpr * N
S <- N - I - ((rho*I*r)/(eta*(1-rho)))
P <- N - I - S

# kappa
I_ambient <- as.vector(I %*% TaR)
kappa <- (I_ambient / H) * c

surv <- peip
Y <- Z / surv
M <- (Z*(g + (f*q*p*kappa))) / (f*q*kappa*p*surv) # M from kappa
lambda <- g*M

M[242] <- 0
Y[242] <- 0
Z[242] <- 0
lambda[242] <- 0
reservoir_eir <- h_patch[242]/b

n_patch <- nrow(TaR)

# --------------------------------------------------------------------------------
#   model state
# --------------------------------------------------------------------------------

# initial conditions
human_state_det <- matrix(data = 0, nrow = n_patch, ncol = 3, dimnames = list(NULL, c('S', 'I', 'P')))
human_state_det[, 'S'] <- S
human_state_det[, 'I'] <- I
human_state_det[, 'P'] <- P

# sample the human state stochastically to get the right total pop
human_state_sto <- human_state_det
for (i in 1:nrow(human_state_sto)) {
  if (N[i] == 0) {
    next
  }
  human_state_sto[i, ] <- as.vector(rmultinom(n = 1, size = N[i], prob = human_state_det[i, ]))
}

MYZ_sto <- matrix(data = 0, nrow = n_patch, ncol = 3, dimnames = list(NULL, c('M', 'Y', 'Z')))
MYZ_sto[, 'M'] <- as.vector(rmultinom(n = 1, size = round(sum(M)), prob = M))
Y_prob <- Y/M
Y_prob[242] <- 0
MYZ_sto[, 'Y'] <- rbinom(n = nrow(MYZ_sto), size = MYZ_sto[, 'M'], prob = Y_prob)
Z_prob <- Z/Y
Z_prob[242] <- 0
MYZ_sto[, 'Z'] <- rbinom(n = nrow(MYZ_sto), size = MYZ_sto[, 'Y'], prob = Z_prob)

tmax <- 365 * 10

mod_det <- make_MicroMoB(tmax = tmax, p = n_patch)
setup_aqua_trace(model = mod_det, lambda = lambda, stochastic = FALSE)
setup_mosquito_RM(model = mod_det, stochastic = FALSE, f = f, q = q, eip = eip, p = p, psi = diag(n_patch), nu = 25, M = M, Y = Y, Z = Z)
setup_humans_SIP(model = mod_det, stochastic = FALSE, theta = TaR, SIP = human_state_det, b = b, c = c, r = r, rho = rho, eta = eta)

mod_sto <- make_MicroMoB(tmax = tmax, p = n_patch)
setup_aqua_trace(model = mod_sto, lambda = lambda, stochastic = TRUE)
setup_mosquito_RM(model = mod_sto, stochastic = TRUE, f = f, q = q, eip = eip, p = p, psi = diag(n_patch), nu = 25, M = MYZ_sto[, 'M'], Y = MYZ_sto[, 'Y'], Z = MYZ_sto[, 'Z'])
setup_humans_SIP(model = mod_sto, stochastic = TRUE, theta = TaR, SIP = human_state_sto, b = b, c = c, r = r, rho = rho, eta = eta)


# compute deterministic bloodmeal
n <- mod_det$global$n
p <- mod_det$global$p
empty_vec <- rep(0, p)

bm_det <- list()

bm_det$H <- compute_H(mod_det)
bm_det$x <- compute_x(mod_det)
bm_det$wf <- compute_wf(mod_det)
bm_det$Psi <- compute_Psi(mod_det)
bm_det$W <- as.vector(t(bm_det$Psi) %*% (bm_det$wf * bm_det$H))

bm_det$beta <- diag(bm_det$wf, nrow = n, ncol = n) %*% bm_det$Psi %*% diag(1/bm_det$W, nrow = p, ncol = p)

bm_det$f <- compute_f(mod_det, B = empty_vec)
bm_det$q <- compute_q(mod_det, W = W, Wd = empty_vec, B = empty_vec)

bm_det$Z <- compute_Z(mod_det)

bm_det$EIR <- as.vector(bm_det$beta %*% (bm_det$f*bm_det$q*bm_det$Z))
bm_det$kappa <- as.vector(t(bm_det$beta) %*% (bm_det$x*bm_det$H))


# compute stochastic bloodmeal
bm_sto <- list()

bm_sto$H <- compute_H(mod_sto)
bm_sto$x <- compute_x(mod_sto)
bm_sto$wf <- compute_wf(mod_sto)
bm_sto$Psi <- compute_Psi(mod_sto)
bm_sto$W <- as.vector(t(bm_sto$Psi) %*% (bm_sto$wf * bm_sto$H))

bm_sto$beta <- diag(bm_sto$wf, nrow = n, ncol = n) %*% bm_sto$Psi %*% diag(1/bm_sto$W, nrow = p, ncol = p)

bm_sto$f <- compute_f(mod_sto, B = empty_vec)
bm_sto$q <- compute_q(mod_sto, W = W, Wd = empty_vec, B = empty_vec)

bm_sto$Z <- compute_Z(mod_sto)

bm_sto$EIR <- as.vector(bm_sto$beta %*% (bm_sto$f*bm_sto$q*bm_sto$Z))
bm_sto$kappa <- as.vector(t(bm_sto$beta) %*% (bm_sto$x*bm_sto$H))
















