rm(list=ls());gc()

# Load Libraries ####
library(data.table)
library(Matrix)
library(MASS)
library(here)
library(MicroMoB)
library(progress)
library(ggplot2)

# transmission parameters
b = 0.55
c = 0.15

# human infection parameters
FeverPf = 0.1116336
TreatPf = 0.602
r = 1/200 # rate at which people become cured
eta = 1/32 # rate at which prophylaxis wears off
rho = FeverPf*TreatPf # probability of clearing infection through treatment cascade


# 2 patch
pfpr <- c(0.05, 0.2)
N <- c(5e4, 2e4)

odds <- r/(1-rho)*pfpr/(1-(1+rho*r/eta/(1-rho))*pfpr)

TaR <- matrix(data = c(
  0.95, 0.05,
  0.15, 0.85
), nrow = 2, ncol = 2, byrow = TRUE)

# Force of Infection (FOI) in each patch (not on each strata)
h_patch <- as.vector(MASS::ginv(TaR) %*% odds)
h_strata <- as.vector(TaR %*% h_patch)
eir_patch <- h_patch/b
eir_strata <- as.vector(TaR %*% eir_patch)

# human pop
I <- pfpr * N
S <- N - I - ((rho*I*r)/(eta*(1-rho)))
P <- N - I - S


# --------------------------------------------------------------------------------
#   deterministic simulation
# --------------------------------------------------------------------------------

# initial conditions
n_patch <- nrow(TaR)
human_state <- matrix(data = 0, nrow = n_patch, ncol = 3, dimnames = list(NULL, c('S', 'I', 'P')))
human_state[, 'S'] <- S
human_state[, 'I'] <- I
human_state[, 'P'] <- P

tmax <- 365 * 2
mod <- make_MicroMoB(tmax = tmax, p = n_patch)

setup_humans_SIP(model = mod, stochastic = FALSE, theta = TaR, SIP = human_state, b = b, c = c, r = r, rho = rho, eta = eta)

# human output table
human_out <- data.table::CJ(day = 1:tmax, state = c('S', 'I', 'P'), patch = 1:n_patch, value = NaN)
human_out <- human_out[c('S', 'I', 'P'), on="state"]
data.table::setkey(human_out, day, patch)

pb <- progress_bar$new(total = tmax)

# run it
while (get_tnow(mod) <= tmax) {
  mod$human$EIR <- eir_strata
  step_humans(model = mod)
  
  human_out[day==get_tnow(mod) & state=='S', value := mod$human$SIP[, 'S']]
  human_out[day==get_tnow(mod) & state=='I', value := mod$human$SIP[, 'I']]
  human_out[day==get_tnow(mod) & state=='P', value := mod$human$SIP[, 'P']]
  
  mod$global$tnow <- mod$global$tnow + 1L
  pb$tick()
}
human_out[, patch := as.factor(patch)]

ggplot(human_out) +
  geom_line(aes(x = day, y = value, color = patch, group = patch)) +
  facet_wrap(. ~ state, scales = 'free')



# --------------------------------------------------------------------------------
#   stochastic simulation
# --------------------------------------------------------------------------------

tmax <- 365 * 10

parout <- parallel::mclapply(X = 1:20, FUN = function(runid) {
  mod <- make_MicroMoB(tmax = tmax, p = n_patch)
  
  setup_humans_SIP(model = mod, stochastic = TRUE, theta = TaR, SIP = round(human_state), b = b, c = c, r = r, rho = rho, eta = eta)
  
  # human output table
  human_out <- data.table::CJ(day = 1:tmax, state = c('S', 'I', 'P'), patch = 1:n_patch, value = NaN)
  human_out <- human_out[c('S', 'I', 'P'), on="state"]
  data.table::setkey(human_out, day, patch)
  
  # run it
  while (get_tnow(mod) <= tmax) {
    mod$human$EIR <- eir_strata
    step_humans(model = mod)
    
    human_out[day==get_tnow(mod) & state=='S', value := mod$human$SIP[, 'S']]
    human_out[day==get_tnow(mod) & state=='I', value := mod$human$SIP[, 'I']]
    human_out[day==get_tnow(mod) & state=='P', value := mod$human$SIP[, 'P']]
    
    mod$global$tnow <- mod$global$tnow + 1L
  }
  human_out[, patch := as.factor(patch)]
  human_out[, run := as.integer(runid)]
  return(human_out)
}, mc.cores = 10)

parout <- data.table::rbindlist(parout)

sto_mean <- parout[, .('mean' = mean(value)), by = .(state, day, patch)]


ggplot(parout) +
  geom_line(aes(x = day, y = value, color = patch, group = patch), alpha = 0.1) +
  geom_line(data = sto_mean, aes(x = day, y = mean, color = patch)) +
  facet_wrap(. ~ state, scales = 'free')
