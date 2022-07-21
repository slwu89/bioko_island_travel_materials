rm(list = ls()); gc()
library(deSolve)
library(expm)
library(MASS)
library(data.table)
library(ggplot2)

# --------------------------------------------------------------------------------
# check forward solutions (start with Lambda)
# --------------------------------------------------------------------------------

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

M_ix <- 1:4
Y_ix <- 5:8
Z_ix <- 9:12

spat_ode <- function(t, y, pars) {
  M <- y[M_ix]
  Y <- y[Y_ix]
  Z <- y[Z_ix]
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

tmax <- spat_out[, max(time)]

# check M
Omega_inv <- MASS::ginv(Omega)
M_analytic <- as.vector(Omega_inv %*% Lambda)
M_simulation <- spat_out[time == tmax, value[M_ix]]

rbind(M_analytic, M_simulation)

# check Y
fqk_Omega_inv <- MASS::ginv(diag(f*q*kappa) + Omega)
Y_analytic <- as.vector(fqk_Omega_inv %*% (f*q*kappa*M_analytic))
Y_simulation <- spat_out[time == tmax, value[Y_ix]]

rbind(Y_analytic, Y_simulation)

# check Z
Z_analytic <- as.vector((Omega_inv %*% OmegaEIP) %*% (matrix(f*q*kappa) * ((M_analytic) - (fqk_Omega_inv %*% (f*q*kappa*(M_analytic))))))
Z_simulation <- spat_out[time == tmax, value[Z_ix]]

rbind(Z_analytic, Z_simulation)


# --------------------------------------------------------------------------------
# check backwards solutions (start with Z)
# --------------------------------------------------------------------------------

# check M-Y
Z <- Z_simulation
MY_analytic <- as.vector((MASS::ginv(OmegaEIP) %*% Omega %*% Z) / matrix(f*q*kappa))
MY_simulation <- spat_out[time == tmax, value[M_ix] - value[Y_ix]]

rbind(MY_analytic, MY_simulation)

# check Y
Y_analytic <- as.vector(Omega_inv %*% (matrix(f*q*kappa) * matrix(MY_analytic)))

rbind(Y_analytic, Y_simulation)

# check M
M_analytic <- MY_analytic + Y_analytic

rbind(M_analytic, M_simulation)

# check Lambda
Lambda_analytic <- as.vector(Omega %*% M_analytic)

rbind(Lambda_analytic, Lambda)

# # check M
# Omega_inv <- MASS::ginv(Omega)
# M <- as.vector(Omega_inv %*% Lambda)
# M
# spat_out[time == 365, value[1:4]]
# 
# # check Y
# spat_out[time == 365, value[5:8]]
# as.vector(MASS::ginv(diag(f*q*kappa) + Omega) %*% (diag(f*q*kappa) %*% M))
# as.vector(MASS::ginv(diag(f*q*kappa) + Omega) %*% (f*q*kappa*M))
# 
# # check M-Y
# Z <- spat_out[time == 365, value[9:12]]
# as.vector((MASS::ginv(OmegaEIP) %*% Omega %*% Z) / matrix(f*q*kappa))
# spat_out[time == 365, value[1:4]] - spat_out[time == 365, value[5:8]]
# 
# # check Y other way
# M_Y <- as.vector((MASS::ginv(OmegaEIP) %*% Omega %*% Z) / matrix(f*q*kappa))
# 
# as.vector(Omega_inv %*% (matrix(f*q*kappa) * matrix(M_Y)))
# spat_out[time == 365, value[5:8]]
# 
# Y <- as.vector(Omega_inv %*% MASS::ginv(OmegaEIP) %*% Omega %*% Z)
# 
# # check Z the ugly way
# (Omega_inv %*% OmegaEIP) %*% (matrix(f*q*kappa) * ((Omega_inv %*% Lambda) - (MASS::ginv(diag(f*q*kappa) + Omega) %*% (f*q*kappa*(Omega_inv%*%Lambda)))))
# 
# as.vector((Omega_inv %*% OmegaEIP) %*% (matrix(f*q*kappa) * ((Omega_inv %*% Lambda) - (MASS::ginv(diag(f*q*kappa) + Omega) %*% (f*q*kappa*(Omega_inv%*%Lambda))))))
# spat_out[time == 365, value[9:12]]