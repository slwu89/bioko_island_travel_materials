rm(list=ls());gc()

# Load Libraries ####
library(data.table)
library(Matrix)
library(MASS)
library(here)
library(MicroMoB)
library(progress)
library(ggplot2)
library(parallel)

# for spatial dde
library(deSolve)
library(expm)

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


# --------------------------------------------------------------------------------
#   parameters
# --------------------------------------------------------------------------------

# parameters
n <- 242
p <- 242

# mosquito parameters
f <- 0.3
q <- 0.9
a <- f*q
g <- 1/10
tau <- 11
sigma <- 0

b <- 0.55
c <- 0.15

calK <- diag(p)*0

Omega <- diag(g, p) + (diag(sigma, p) %*% (diag(p) - calK))
Omega_inv <- solve(Omega)
OmegaEIP <- expm::expm(-Omega * tau)
OmegaEIP_inv <- expm::expm(Omega * tau)

# transmission parameters
b <- 0.55
c <- 0.15

# human infection parameters
FeverPf <- 0.1116336
TreatPf <- 0.602
r <- 1/200 # rate at which people become cured
eta <- 1/32 # rate at which prophylaxis wears off
rho <- FeverPf*TreatPf # probability of clearing infection through treatment cascade


# --------------------------------------------------------------------------------
# equilibrium solutions
# --------------------------------------------------------------------------------

# human movement
Psi <- t(TaR)
W <- Psi %*% N
beta <- t(Psi) %*% diag(as.vector(1/W), p, p)

odds <- r/(1-rho)*pfpr/(1-(1+rho*r/eta/(1-rho))*pfpr)

# Force of Infection (FOI) in each patch (not on each strata)
h_patch <- as.vector(solve(TaR) %*% odds)

# EIR in patches
H <- as.vector(N %*% TaR)
Z <- (h_patch * H) / (a * b)

I <- pfpr * N
S <- N - I - ((rho*I*r)/(eta*(1-rho)))
P <- N - I - S

# we expect this to be about the same
as.vector(TaR %*% h_patch) - (eta * P) / (rho * S)

# kappa
I_ambient <- as.vector(I %*% TaR)
kappa <- (I_ambient / H) * c

# once we have equilibrium Z, calculate the rest of the mosquito model at equilibrium

# derived M-Y
MY <- diag(1/as.vector(f*q*kappa), p, p) %*% OmegaEIP_inv %*% Omega %*% Z

# derived Y
Y <- Omega_inv %*% (diag(as.vector(f*q*kappa), p, p) %*% MY)

# derived M
M <- MY + Y

# derived Lambda
Lambda <- Omega %*% M

Z <- as.vector(Z)
Y <- as.vector(Y)
M <- as.vector(M)
Lambda <- as.vector(Lambda)

M[242] <- 0
Y[242] <- 0
Z[242] <- 0
Lambda[242] <- 0
rio_muni_eir <- h_patch[242]/b

# parameter set
params <- list(
  b = b,
  c = c,
  r = r,
  f = f,
  q = q,
  g = g,
  sigma = sigma,
  tau = tau,
  p = p,
  n = n,
  Omega = Omega,
  OmegaEIP = OmegaEIP,
  beta = beta,
  betaT = t(beta),
  Psi = Psi,
  N = N,
  eta = eta,
  rho = rho,
  rio_muni_eir = rio_muni_eir
)

make_index <- function(pars, L = FALSE) {
  n <- pars$n
  p <- pars$p
  if (L) {
    pars$L_ix <- 1:p
    max_ix <- p
  } else {
    pars$L_ix <- integer(0)
    max_ix <- 0
  }
  pars$M_ix <- seq(from = max_ix+1, length.out = p)
  max_ix <- tail(pars$M_ix, 1)
  
  pars$Y_ix <- seq(from = max_ix+1, length.out = p)
  max_ix <- tail(pars$Y_ix, 1)
  
  pars$Z_ix <- seq(from = max_ix+1, length.out = p)
  max_ix <- tail(pars$Z_ix, 1)
  
  pars$X_ix <- seq(from = max_ix+1, length.out = n)
  max_ix <- tail(pars$X_ix, 1)
  
  pars$P_ix <- seq(from = max_ix+1, length.out = n)
  
  return(pars)
}

params <- make_index(params)


# --------------------------------------------------------------------------------
#   spatial DDEs
# --------------------------------------------------------------------------------

# compute kappa at time t
compute_kappa <- function(y, pars) {
  pars$betaT %*% (y[pars$X_ix] * pars$c)
}

# compute kappa at time t-tau
compute_kappa_tau <- function(t, y, pars) {
  if (t < pars$tau) {
    X_tau <- pars$Y0[pars$X_ix]
  } else {
    X_tau <- deSolve::lagvalue(t = t - pars$tau, nr = pars$X_ix)
  }
  pars$betaT %*% (X_tau * pars$c)
}

# basic spatial differential equations for RM mosquito model
mosquito_dt <- function(t, y, pars, Lambda, kappa, kappa_tau) {
  M <- y[pars$M_ix]
  Y <- y[pars$Y_ix]
  Z <- y[pars$Z_ix]
  if (t < pars$tau) {
    M_tau <- pars$Y0[pars$M_ix]
    Y_tau <- pars$Y0[pars$Y_ix]
  } else {
    M_tau <- deSolve::lagvalue(t = t - tau, nr = pars$M_ix)
    Y_tau <- deSolve::lagvalue(t = t - tau, nr = pars$Y_ix)
  }
  dMdt <- Lambda - (pars$Omega %*% M)
  dYdt <- diag(as.vector(pars$f*pars$q*kappa)) %*% (M - Y) - (pars$Omega %*% Y)
  dZdt <- OmegaEIP %*% diag(as.vector(pars$f*pars$q*kappa_tau)) %*% (M_tau - Y_tau) - (pars$Omega %*% Z)
  return(c(dMdt, dYdt, dZdt))
}


# lambda at time t for trace-based aquatic mosquito model
compute_Lambda.null <- function(y, pars) {
  return(pars$Lambda)
}

# null differential equations for trace-based aquatic mosquito model
aquatic_dt.null <- function(t, y, pars) {
  return(numeric(0))
}

# compute EIR
compute_EIR <- function(y, pars) {
  Z <- y[pars$Z_ix]
  pars$beta %*% diag(pars$f * pars$q, nrow = pars$p) %*% Z
}

# basic spatial SIP (Susceptible-Infectious-Protected) model
human_dt <- function(t, y, pars, EIR) {
  X <- y[pars$X_ix]
  P <- y[pars$P_ix]
  dXdt <- (diag((1-pars$rho)*pars$b*EIR, nrow = pars$n, ncol = pars$n) %*% (pars$N - X - P)) - (pars$r * X)
  dYdt <- (diag(pars$rho*pars$b*EIR, nrow = pars$n, ncol = pars$n) %*% (pars$N - X - P)) - (pars$eta * P)
  return(c(dXdt, dYdt))
}


# complete spatial DDE model
spatialdynamics_dt <- function(t, y, pars) {
  
  # bloodfeeding
  EIR <- compute_EIR(y, pars)
  EIR <- as.vector(EIR)
  EIR <- EIR + (pars$Psi[242, ] * pars$rio_muni_eir)
  
  kappa <- compute_kappa(y, pars)
  kappa_tau <- compute_kappa_tau(t, y, pars)
  
  # aquatic stage
  dL <- aquatic_dt.null(t, y, pars)
  Lambda <- compute_Lambda.null(y, pars)
  
  # adult stage
  dMYZ <- mosquito_dt(t, y, pars, Lambda, kappa, kappa_tau)
  
  # human
  dX <- human_dt(t, y, pars, EIR)
  
  return(list(c(dL, dMYZ, dX)))
  
}


# --------------------------------------------------------------------------------
#   deterministic simulation
# --------------------------------------------------------------------------------


Y0 <- rep(NaN, max(params$P_ix))
Y0[params$M_ix] <- M
Y0[params$Y_ix] <- Y
Y0[params$Z_ix] <- Z
Y0[params$X_ix] <- I
Y0[params$P_ix] <- P

params$Y0 <- Y0
params$Lambda <- Lambda

system.time(
  out <- deSolve::dede(y = params$Y0, times = 0:(365*2), func = spatialdynamics_dt, parms = params)
)

# check we're at equilibrium
rbind(as.numeric(out[nrow(out), params$X_ix+1]), I)
