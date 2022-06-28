rm(list=ls());gc()

# --------------------------------------------------------------------------------
# figuring out how to set up the equilibrium model in 3 patch case
# --------------------------------------------------------------------------------

# TaR <- matrix(
#   data = c(
#     0.9, 0.025, 0.075,
#     0.1, 0.85, 0.05,
#     0.01, 0.09, 0.9
#   ),
#   nrow = 3, ncol = 3, byrow = TRUE
# )
# 
# # pfpr amongst residents
# pfpr <- c(0.4, 0.25, 0.1)
# 
# # total human pop
# H <- c(100, 200, 300)
# 
# 
# 
# 
# 
# # See Ruktanonchai et al. (2016) and Supplementary Information for derivation
# odds.vector <- r/(1-rho)*pfpr/(1-(1+rho*r/eta/(1-rho))*pfpr)
# 
# # Force of Infection (FOI) or "happenings rate" h
# h.FOI <- MASS::ginv(TaR.matrix) %*% odds.vector
# h.FOI[which(h.FOI < 0)] <- 0
# h.FOI <- as.vector(h.FOI)
# 
# # Calculate Mosquito Parameters ####
# # we begin by calculating the fraction of infected mosquitoes, the sporozoite rate, from the Ross-Macdonald equations + Kappa
# # Total population, including visitors
# H.visitors <- t(pop.data[year == 2018, pop] %*% TaR.matrix)
# # Sick population, including visitors
# X.visitors <- t((pop.data[year == 2018, pop] * pfpr)  %*% TaR.matrix)
# 
# # kappa: net infectiousness of human pop at each place
# kappa <- X.visitors/H.visitors * c
# kappa <- as.vector(kappa)
# 
# H.visitors <- as.vector(H.visitors)
# n_patch <- length(kappa)


# --------------------------------------------------------------------------------
# test SIP equilibrium for non spatial model
# --------------------------------------------------------------------------------

library(deSolve)

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

pfpr <- 0.2
N <- 500

# eq equations
I <- pfpr * N
S <- N - I - ((rho*I*r)/(eta*(1-rho)))
P <- N - I - S
h <- (eta * P) / (rho * S)

sip_ode <- function(t, y, params) {
  eta <- params$eta
  r <- params$r
  h <- params$h
  rho <- params$rho
  dxdt <- rep(0, 3)
  S <- y[1]
  I <- y[2]
  P <- y[3]
  dxdt[1] <- -h*S + r*I + eta*P
  dxdt[2] <- h*(1-rho)*S - r*I
  dxdt[3] <- h*rho*S - eta*P
  return(list(dxdt))
}

params <- list(r=r,h=h,rho=rho,eta=eta)
y <- c(S,I,P)

out <- ode(y = y, times = 0:100, func = sip_ode, parms = params, method = "ode23")


# --------------------------------------------------------------------------------
# testing FOI in space
# --------------------------------------------------------------------------------

TaR <- matrix(
  data = c(
    0.9, 0.025, 0.075,
    0.1, 0.85, 0.05,
    0.01, 0.09, 0.9
  ),
  nrow = 3, ncol = 3, byrow = TRUE
)

N <- c(100, 200, 300) #resident pop
pfpr <- c(0.1, 0.3, 0.25)

# EIR <- c(0.5, 2, 1)
# 
# h_manual <- rep(0, 3)
# for (i in 1:3) {
#   for (j in 1:3) {
#     h_manual[i] <- h_manual[i] + (TaR[i,j] * EIR[j])
#   }
# }
# 
# h_mat <- TaR %*% EIR
# 
# 
# H_man <- rep(0, 3)
# for (j in 1:3) {
#   for (i in 1:3) {
#     H_man[j] <- H_man[j] + (TaR[i,j] * N[i])
#   }
# }
# 
# H_mat <- N %*% TaR

# let's do this 3 patch setup and try to solve the EIR required to give equilibrium PfPR
I <- pfpr * N
S <- N - I - ((rho*I*r)/(eta*(1-rho)))
P <- N - I - S
h <- (eta * P) / (rho * S) # the FOI experienced by each population during their movement

EIR <- MASS::ginv(TaR) %*% h/b # the EIR produced by mosquitoes at each patch
EIR <- as.vector(EIR)

# check, should be very small
abs(h - as.vector(TaR %*% EIR)*b)

# get Z
H <- N %*% TaR
H <- as.vector(H)
Z <- (EIR * H) / a # Z calc from ambient pop, not census pop

# kappa
c <- 0.15
# kappa <- (pfpr*c) %*% TaR
# kappa <- as.vector(kappa)

# kappa a diff way (correct)
# ambient pfpr
kappa <- as.vector((I %*% TaR) / (N %*% TaR)) * c

# rest of the mosy
surv <- peip
Y <- Z / surv
M <- (Z*(g + (f*q*p*kappa))) / (f*q*kappa*p*surv) # M from kappa
lambda <- g*M


# --------------------------------------------------------------------------------
# bloodmeal calcs
# --------------------------------------------------------------------------------

bm <- list()

# human quantities
bm$H <- N
bm$x <- pfpr * c
bm$wf <- rep(1,3)
bm$Psi <- TaR
bm$W <- as.vector(t(bm$Psi) %*% (bm$wf * bm$H))

# biting distribution matrix (n x p)
bm$beta <- diag(bm$wf, nrow = 3, ncol = 3) %*% bm$Psi %*% diag(1/bm$W, nrow = 3, ncol = 3)

# density of infective mosquitoes
bm$Z <- Z

# EIR experienced by strata
bm$EIR <- as.vector(bm$beta %*% (a*bm$Z))

# kappa on mosquitoes from ambient pop
bm$kappa <- as.vector(t(bm$beta) %*% (bm$x*bm$H))






# --------------------------------------------------------------------------------
#   deterministic simulation
# --------------------------------------------------------------------------------

library(MicroMoB)
library(progress)
library(ggplot2)

# initial conditions
n_patch <- 3
human_state <- matrix(data = 0, nrow = n_patch, ncol = 3, dimnames = list(NULL, c('S', 'I', 'P')))
human_state[, 'S'] <- S
human_state[, 'I'] <- I
human_state[, 'P'] <- P

tmax <- 365 * 2
mod <- make_MicroMoB(tmax = tmax, p = n_patch)

setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE)
setup_mosquito_RM(model = mod, stochastic = FALSE, f = f, q = q, eip = eip, p = p, psi = diag(3), nu = 25, M = M, Y = Y, Z = Z)
setup_humans_SIP(model = mod, stochastic = FALSE, theta = TaR, SIP = human_state, b = b, c = c, r = r, rho = rho, eta = eta)
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

