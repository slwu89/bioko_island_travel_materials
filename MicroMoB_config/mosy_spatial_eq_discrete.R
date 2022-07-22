rm(list=ls()); gc()
library(deSolve)
library(expm)
library(MASS)
library(data.table)
library(ggplot2)
library(MicroMoB)

p <- 1-(1/20)
g <- 1 - p

f <- 0.3
q <- 0.9

n_patch <- 4

tau <- 11

sigma <- 1/20 # daily prob of leaving

# for sampling K
set.seed(43583491L)
K <- matrix(0, n_patch, n_patch)
K[upper.tri(K)] <- rexp(sum(1:(n_patch-1)))
K[lower.tri(K)] <- rexp(sum(1:(n_patch-1)))

Psi <- K
Psi <- Psi/(rowSums(Psi)/sigma)
diag(Psi) <- 1 - sigma

K_analytic <- K/rowSums(K)

kappa <- c(0.1, 0.075, 0.025, 0.05)

Lambda <- c(5, 10, 8, 6)



tmax <- 500

M <- rep(0, 4)
Y <- rep(0, 4)
Z <- rep(0, 4)

mod <- make_MicroMoB(tmax = tmax, p = 4)
setup_mosquito_RM(mod, stochastic = FALSE, f = f, q = q, eip = tau, p = p, psi = Psi, M = M, Y = Y, Z = Z)
setup_aqua_trace(model = mod, lambda = Lambda, stochastic = FALSE)

det_out <- data.table::CJ(day = 1:tmax, state = c('M', 'Y', 'Z'), patch = 1:4, value = NaN)
det_out <- det_out[c('M', 'Y', 'Z'), on="state"]
data.table::setkey(det_out, day)

# run it
while(get_tnow(mod) <= tmax) {
  
  mod$mosquito$kappa <- kappa
  step_aqua(model = mod)
  step_mosquitoes(model = mod)
  
  det_out[day == get_tnow(mod) & state == 'M', value := mod$mosquito$M]
  det_out[day == get_tnow(mod) & state == 'Y', value := mod$mosquito$Y]
  det_out[day == get_tnow(mod) & state == 'Z', value := mod$mosquito$Z]
  
  mod$global$tnow <- mod$global$tnow + 1L
}

ggplot(det_out) +
  geom_line(aes(x=day,y=value,color=state)) +
  facet_grid(state ~ patch, scales = 'free')




Omega <- diag(g, n_patch) + (diag(sigma, n_patch) %*% (diag(n_patch) - K_analytic))
OmegaEIP <- expm::expm(-Omega * tau)
Omega_inv <- solve(Omega)

# check M
M_analytic <- as.vector(Omega_inv %*% Lambda)
M_simulation <- det_out[day == tmax & state == 'M', value]

rbind(M_analytic, M_simulation)

# M discrete time
M_analytic <- t(Lambda) %*% solve(diag(4) - (Psi %*% diag(p, 4)))

# check Y
# fqk_Omega_inv <- MASS::ginv(diag(f*q*kappa) + Omega)
fqk_Omega_inv <- solve(diag(f*q*kappa) + Omega)
Y_analytic <- as.vector(fqk_Omega_inv %*% (f*q*kappa*M_analytic))
Y_simulation <- det_out[day == tmax & state == 'Y', value]

rbind(Y_analytic, Y_simulation)

# Y discrete time

Y_analytic <- (M_analytic %*% diag(f*q*kappa) %*% Psi %*% diag(p,4)) %*% solve(diag(4) - (Psi %*% diag(p,4)) + (diag(f*q*kappa) %*% Psi %*% diag(p,4)))


# check Z
Z_analytic <- as.vector((Omega_inv %*% OmegaEIP) %*% (matrix(f*q*kappa) * ((M_analytic) - (fqk_Omega_inv %*% (f*q*kappa*(M_analytic))))))
Z_simulation <- det_out[day == tmax & state == 'Z', value]

rbind(Z_analytic, Z_simulation)

# Z discrete time

Z_analytic <- ((M_analytic - Y_analytic) %*% diag(f*q*kappa) %*% (Psi%^%(tau+1)) %*% (diag(p,4)%^%(tau+1))) %*% solve(diag(4) - (Psi %*% diag(p,4)))

# M-Y discrete time (backwards)
M_Y_analytic <- Z_analytic %*% (diag(4) - (Psi %*% diag(p,4))) %*% solve(diag(f*q*kappa) %*% (Psi%^%(tau+1)) %*% (diag(p,4)%^%(tau+1)))

# Y discrete time (backwards)
M_Y_analytic %*% diag(f*q*kappa) %*% Psi %*% diag(p,4) %*% solve(diag(4) - (Psi %*% diag(p,4)))

# Lambda discrete time (backwards)
M_analytic %*% (diag(4) - (Psi %*% diag(p,4)))
