rm(list=ls());gc()

# --------------------------------------------------------------------------------
# figuring out how to set up the equilibrium model in 3 patch case
# --------------------------------------------------------------------------------

TaR <- matrix(
  data = c(
    0.9, 0.025, 0.075,
    0.1, 0.85, 0.05,
    0.01, 0.09, 0.9
  ),
  nrow = 3, ncol = 3, byrow = TRUE
)

# pfpr amongst residents
pfpr <- c(0.4, 0.25, 0.1)

# total human pop
H <- c(100, 200, 300)





# See Ruktanonchai et al. (2016) and Supplementary Information for derivation
odds.vector <- r/(1-rho)*pfpr/(1-(1+rho*r/eta/(1-rho))*pfpr)

# Force of Infection (FOI) or "happenings rate" h
h.FOI <- MASS::ginv(TaR.matrix) %*% odds.vector
h.FOI[which(h.FOI < 0)] <- 0
h.FOI <- as.vector(h.FOI)

# Calculate Mosquito Parameters ####
# we begin by calculating the fraction of infected mosquitoes, the sporozoite rate, from the Ross-Macdonald equations + Kappa
# Total population, including visitors
H.visitors <- t(pop.data[year == 2018, pop] %*% TaR.matrix)
# Sick population, including visitors
X.visitors <- t((pop.data[year == 2018, pop] * pfpr)  %*% TaR.matrix)

# kappa: net infectiousness of human pop at each place
kappa <- X.visitors/H.visitors * c
kappa <- as.vector(kappa)

H.visitors <- as.vector(H.visitors)
n_patch <- length(kappa)



# test SIP equilibrium for non spatial model
library(deSolve)

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
