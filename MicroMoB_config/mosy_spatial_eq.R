library(deSolve)
library(expm)
library(MASS)
library(data.table)
library(ggplot2)

f <- 0.3
q <- 0.9

n_patch <- 4
g <- 1/20
sigma <- 1/10

tau <- 11

K <- matrix(0, n_patch, n_patch)
K[upper.tri(K)] <- rexp(sum(1:(n_patch-1)))
K[lower.tri(K)] <- rexp(sum(1:(n_patch-1)))
K <- K/rowSums(K)

Omega <- diag(g, n_patch) + (diag(sigma, n_patch) %*% (diag(n_patch) - K))
OmegaEIP <- expm::expm(-Omega * tau)

kappa <- c(0.1, 0.075, 0.025, 0.05)

Lambda <- c(5, 10, 8, 6)

Y0 <- rep(0, n_patch*3)

spat_ode <- function(t, y, pars) {
  M <- y[1:4]
  Y <- y[5:8]
  Z <- y[9:12]
  if (t < tau) {
    M_tau <- Y0[1:4]
    Y_tau <- Y0[5:8]
  } else {
    MYZ_tau <- lagvalue(t - tau)
    M_tau <- MYZ_tau[1:4]
    Y_tau <- MYZ_tau[5:8]
  }
  dMdt <- Lambda - (Omega %*% M)
  dYdt <- f*q*kappa*(M - Y) - (Omega %*% Y)
  dZdt <- (OmegaEIP %*% (f*q*kappa*(M_tau - Y_tau))) - (Omega %*% Z)
  return(list(c(dMdt, dYdt, dZdt)))
}

spat_out <- dede(y = Y0, times = 0:365, func = spat_ode)
spat_out <- as.data.frame(spat_out)
colnames(spat_out) <- c("time", unlist(lapply(c("M","Y","Z"), function(x) paste0(x, 1:4))))

spat_out <- melt(as.data.table(spat_out), id.vars = 'time')

ggplot(spat_out) +
  geom_line(aes(x=time,y=value,color=variable)) +
  facet_wrap(. ~ variable, scales = 'free')

# check results in SI
Omega_inv <- MASS::ginv(Omega)
Omega_inv %*% Lambda
spat_out[time == 365, value[1:4]]

# solve the inverse, given Z solve for Lambda to give that, and the other state variables
Z <- c(50, 20, 100, 80)
OmegaEIP_inv <- MASS::ginv(OmegaEIP)

OmegaEIP_inv %*% Omega %*% Z

as.vector((OmegaEIP_inv %*% (Omega %*% Z)) / (f*q*kappa))

as.vector(OmegaEIP_inv %*% (Omega %*% Z))




