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

# for sampling K
set.seed(43583491L)
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
    M_tau <- Y0[M_ix]
    Y_tau <- Y0[Y_ix]
  } else {
    MYZ_tau <- lagvalue(t - tau)
    M_tau <- MYZ_tau[M_ix]
    Y_tau <- MYZ_tau[Y_ix]
  }
  dMdt <- Lambda - (Omega %*% M)
  dYdt <- diag(f*q*kappa) %*% (M - Y) - (Omega %*% Y)
  dZdt <- OmegaEIP %*% diag(f*q*kappa) %*% (M_tau - Y_tau) - (Omega %*% Z)
  return(list(c(dMdt, dYdt, dZdt)))
}

spat_out <- dede(y = Y0, times = 0:500, func = spat_ode)
spat_out <- as.data.frame(spat_out)
colnames(spat_out) <- c("time", unlist(lapply(c("M","Y","Z"), function(x) paste0(x, 1:4))))

spat_out <- melt(as.data.table(spat_out), id.vars = 'time')

ggplot(spat_out) +
  geom_line(aes(x=time,y=value,color=variable)) +
  facet_wrap(. ~ variable, scales = 'free')

tmax <- spat_out[, max(time)]

# useful matrices
# OmegaEIP_inv <- qr.solve(OmegaEIP)
OmegaEIP_inv <- expm::expm(Omega * tau) # properties of matrix exp
Omega_inv <- solve(Omega)

# check M
M_analytic <- as.vector(Omega_inv %*% Lambda)
M_simulation <- spat_out[time == tmax, value[M_ix]]

rbind(M_analytic, M_simulation)

# check Y
# fqk_Omega_inv <- MASS::ginv(diag(f*q*kappa) + Omega)
fqk_Omega_inv <- solve(diag(f*q*kappa) + Omega)
# Y_analytic <- as.vector(fqk_Omega_inv %*% (f*q*kappa*M_analytic))
Y_simulation <- spat_out[time == tmax, value[Y_ix]]

Y_analytic <- as.vector(fqk_Omega_inv %*% diag(f*q*kappa) %*% matrix(M_analytic))


rbind(Y_analytic, Y_simulation)

# check Z
Z_analytic <- as.vector((Omega_inv %*% OmegaEIP) %*% (diag(f*q*kappa) %*% ((M_analytic) - (fqk_Omega_inv %*% (diag(f*q*kappa) %*% matrix(M_analytic))))))
Z_simulation <- spat_out[time == tmax, value[Z_ix]]

rbind(Z_analytic, Z_simulation)


# --------------------------------------------------------------------------------
# check backwards solutions (start with Z)
# --------------------------------------------------------------------------------

# check M-Y
Z <- Z_simulation
# MY_analytic <- as.vector((OmegaEIP_inv %*% Omega %*% Z) / matrix(f*q*kappa))
MY_analytic <- as.vector(diag(1/(f*q*kappa)) %*% OmegaEIP_inv %*% Omega %*% Z)
MY_simulation <- spat_out[time == tmax, value[M_ix] - value[Y_ix]]

rbind(MY_analytic, MY_simulation)

# check Y
Y_analytic <- as.vector(Omega_inv %*% (diag(f*q*kappa) %*% matrix(MY_analytic)))

rbind(Y_analytic, Y_simulation)

# check M
M_analytic <- MY_analytic + Y_analytic

rbind(M_analytic, M_simulation)

# check Lambda
Lambda_analytic <- as.vector(Omega %*% M_analytic)

rbind(Lambda_analytic, Lambda)
