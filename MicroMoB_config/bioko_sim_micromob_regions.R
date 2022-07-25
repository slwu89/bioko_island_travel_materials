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


# --------------------------------------------------------------------------------
#   aggregate landscape
# --------------------------------------------------------------------------------

clusters <- list("Malabo" = c("Malabo", "Peri"), "Baney" = "Baney", "Luba" = c("Luba", "Moka", "Ureka"), "Riaba" = "Riaba")

# clusters_ix <- vapply(X = pop.data[year == 2018, ad2], FUN = function(i) {
#   which(vapply(X = clusters, FUN = '%in%', FUN.VALUE = logical(1), x = i), useNames = FALSE)
# }, FUN.VALUE = integer(1), USE.NAMES = FALSE)

clusters_ix <- lapply(X = names(clusters), FUN = function(x) {
  which(unlist(lapply(X = pop.data[year == 2018, ad2], FUN = '==', y = x)))
})

clusters_eta <- lapply(X = 1:4, FUN = function(x) {
  pop.data[year == 2018, pop[clusters_ix[[x]]]] / pop.data[year == 2018, sum(pop[clusters_ix[[x]]])]
})

TaR_cluster <- matrix(data = 0, nrow = 4, ncol = 4, dimnames = list(names(clusters), names(clusters)))

# iterate over AI,AJ
for (I in 1:4) {
  for (J in 1:4) {
   # iterate over i and j
    p_IJ <- 0
    for (i in 1:length(clusters_ix[[I]])) {
      for (j in 1:length(clusters_ix[[J]])) {
        p_IJ <- p_IJ + (TaR[clusters_ix[[I]][i], clusters_ix[[J]][j]]  * clusters_eta[[I]][i])
      }
    }
    TaR_cluster[I, J] <- p_IJ
  }
}

TaR_cluster <- TaR_cluster / rowS


# eg_pop <- 1308975
# rio_muni_pop <- eg_pop - sum(N)
# N[242] <- rio_muni_pop

# pfpr <- pfpr.data$draw.mean
# N <- pop.data[year == 2018]$pop
# TaR <- TaR[1:241, ]

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
# equilibrium solutions
# --------------------------------------------------------------------------------

# hack to fix problematically low pfpr values
pfpr[c(71,72,109,215,216)] <- pfpr[c(71,72,109,215,216)] * 1.5
N <- N*1e4
# pfpr[c(71,72,109,215,216)] <- pfpr[c(71,72,109,215,216)] + 0.05
# pfpr[c(71,72,109,215,216)] <- 0

odds <- r/(1-rho)*pfpr/(1-(1+rho*r/eta/(1-rho))*pfpr)

# Force of Infection (FOI) in each patch (not on each strata)
h_patch <- as.vector(MASS::ginv(TaR) %*% odds)

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

surv <- peip
Y <- Z / surv
M <- (Z*(g + (f*q*p*kappa))) / (f*q*kappa*p*surv) # M from kappa
lambda <- g*M

M[242] <- 0
Y[242] <- 0
Z[242] <- 0
lambda[242] <- 0
rio_muni_eir <- h_patch[242]/b


# --------------------------------------------------------------------------------
#   deterministic simulation
# --------------------------------------------------------------------------------

# initial conditions
n_patch <- nrow(TaR)
human_state <- matrix(data = 0, nrow = n_patch, ncol = 3, dimnames = list(NULL, c('S', 'I', 'P')))
human_state[, 'S'] <- S
human_state[, 'I'] <- I
human_state[, 'P'] <- P
human_state[242, ] <- 0

tmax <- 365 * 10
mod <- make_MicroMoB(tmax = tmax, p = n_patch)

setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE)
setup_mosquito_RM(model = mod, stochastic = FALSE, f = f, q = q, eip = eip, p = p, psi = diag(n_patch), nu = 25, M = M, Y = Y, Z = Z)
setup_humans_SIP(model = mod, stochastic = FALSE, theta = TaR, SIP = human_state, b = b, c = c, r = r, rho = rho, eta = eta)

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
  # control
  compute_bloodmeal_simple(model = mod)
  mod$human$EIR <- mod$human$EIR + (TaR[,242] * rio_muni_eir)
  
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

human_out <- human_out[, .(value = sum(value)), by = .(state, day)]
mosy_out <- mosy_out[, .(value = sum(value)), by = .(state, day)]

human_out[, species := "human"]
mosy_out[, species := "mosquito"]
det_out <- rbind(human_out, mosy_out)


# --------------------------------------------------------------------------------
#   stochastic simulation
# --------------------------------------------------------------------------------

parout <- parallel::mclapply(X = 1:40, FUN = function(runid) {
  
  mod <- make_MicroMoB(tmax = tmax, p = n_patch)
  
  human_state_sto <- human_state
  for (i in 1:nrow(human_state_sto)) {
    if (N[i] == 0) {
      next
    }
    human_state_sto[i, ] <- as.vector(rmultinom(n = 1, size = N[i], prob = human_state[i, ]))
  }
  
  M <- round(M)
  Y <- round(Y)
  Z <- round(Z)
  
  setup_aqua_trace(model = mod, lambda = lambda, stochastic = TRUE)
  setup_mosquito_RM(model = mod, stochastic = TRUE, f = f, q = q, eip = eip, p = p, psi = diag(n_patch), nu = 25, M = M, Y = Y, Z = Z)
  setup_humans_SIP(model = mod, stochastic = TRUE, theta = TaR, SIP = human_state_sto, b = b, c = c, r = r, rho = rho, eta = eta)
  
  # human output table
  human_out <- data.table::CJ(day = 1:tmax, state = c('S', 'I', 'P'), patch = 1:n_patch, value = NaN)
  human_out <- human_out[c('S', 'I', 'P'), on="state"]
  data.table::setkey(human_out, day, patch)
  
  # mosy output table
  mosy_out <- data.table::CJ(day = 1:tmax, state = c('M', 'Y', 'Z'), patch = 1:n_patch, value = NaN)
  mosy_out <- mosy_out[c('M', 'Y', 'Z'), on="state"]
  data.table::setkey(mosy_out, day, patch)
  
  # run it
  while (get_tnow(mod) <= tmax) {
    
    compute_bloodmeal_simple(model = mod)
    mod$human$EIR <- mod$human$EIR + (TaR[,242] * rio_muni_eir)
    
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
  }
  
  human_out <- human_out[, .(value = sum(value)), by = .(state, day)]
  mosy_out <- mosy_out[, .(value = sum(value)), by = .(state, day)]
  
  human_out[, species := "human"]
  mosy_out[, species := "mosquito"]
  
  out <- rbind(mosy_out, human_out)
  out[, run := as.integer(runid)]
  
  return(out)
}, mc.cores = 20)

parout <- data.table::rbindlist(parout)

sto_mean <- parout[, .('mean' = mean(value)), by = .(state, day, species)]

ggplot(parout) +
  geom_line(aes(x = day, y = value, color = species, linetype = state, group = run), alpha = 0.05) +
  geom_line(data = det_out, aes(x = day, y = value, color = species)) +
  geom_line(data = sto_mean, aes(x = day, y = mean, color = species)) +
  facet_wrap(species ~ state, scales = "free") +
  ggtitle("summarized Bioko equilibrium simulation")
