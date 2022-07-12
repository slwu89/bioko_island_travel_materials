rm(list=ls());gc()
library(MicroMoB)
library(progress)
library(ggplot2)
library(deSolve)


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

# pfpr <- 0.2
# N <- 500
# 
# # eq equations for ODE
# I <- pfpr * N
# S <- N - I - ((rho*I*r)/(eta*(1-rho)))
# P <- N - I - S
# h <- (eta * P) / (rho * S)
# 
# sip_ode <- function(t, y, params) {
#   eta <- params$eta
#   r <- params$r
#   h <- params$h
#   rho <- params$rho
#   dxdt <- rep(0, 3)
#   S <- y[1]
#   I <- y[2]
#   P <- y[3]
#   dxdt[1] <- -h*S + r*I + eta*P
#   dxdt[2] <- h*(1-rho)*S - r*I
#   dxdt[3] <- h*rho*S - eta*P
#   return(list(dxdt))
# }
# 
# params <- list(r=r,h=h,rho=rho,eta=eta)
# y <- c(S,I,P)
# 
# out <- ode(y = y, times = 0:100, func = sip_ode, parms = params, method = "ode23")
# plot(out)


# --------------------------------------------------------------------------------
# testing FOI in space
# --------------------------------------------------------------------------------

TaR <- matrix(
  data = c(
    0.9, 0.025, 0.07, 0.005,
    0.1, 0.8, 0.05, 0.05,
    0.01, 0.09, 0.875, 0.025,
    0, 0, 0, 1
  ),
  nrow = 4, ncol = 4, byrow = TRUE
)

N <- c(1, 2, 3, 0) * 1e6 #resident pop
pfpr <- c(0.1, 0.3, 0.25, 0.15)

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

M[4] <- 0
Y[4] <- 0
Z[4] <- 0
lambda[4] <- 0
reservoir_eir <- h_patch[4]/b


# --------------------------------------------------------------------------------
#   deterministic simulation
# --------------------------------------------------------------------------------

# initial conditions
n_patch <- nrow(TaR)
human_state <- matrix(data = 0, nrow = n_patch, ncol = 3, dimnames = list(NULL, c('S', 'I', 'P')))
human_state[, 'S'] <- S
human_state[, 'I'] <- I
human_state[, 'P'] <- P

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
  compute_bloodmeal_simple(model = mod)
  mod$human$EIR <- mod$human$EIR + (TaR[,4] * reservoir_eir)

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

ggplot(det_out) +
  geom_line(aes(x = day, y = value, color = species)) +
  facet_wrap(species ~ state, scales = "free")



# --------------------------------------------------------------------------------
#   stochastic simulation
# --------------------------------------------------------------------------------

parout <- parallel::mclapply(X = 1:40, FUN = function(runid) {
  
  mod <- make_MicroMoB(tmax = tmax, p = n_patch)
  
  human_state <- round(human_state)
  M <- round(M)
  Y <- round(Y)
  Z <- round(Z)
  
  setup_aqua_trace(model = mod, lambda = lambda, stochastic = TRUE)
  setup_mosquito_RM(model = mod, stochastic = TRUE, f = f, q = q, eip = eip, p = p, psi = diag(n_patch), nu = 25, M = M, Y = Y, Z = Z)
  setup_humans_SIP(model = mod, stochastic = TRUE, theta = TaR, SIP = human_state, b = b, c = c, r = r, rho = rho, eta = eta)
  
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
    mod$human$EIR <- mod$human$EIR + (TaR[,4] * reservoir_eir)

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
}, mc.cores = 10)

parout <- data.table::rbindlist(parout)

sto_mean <- parout[, .('mean' = mean(value)), by = .(state, day, species)]

ggplot(parout) +
  geom_line(aes(x = day, y = value, color = species, linetype = state, group = run), alpha = 0.05) +
  geom_line(data = det_out, aes(x = day, y = value, color = species)) +
  geom_line(data = sto_mean, aes(x = day, y = mean, color = species)) +
  facet_wrap(species ~ state, scales = "free")

