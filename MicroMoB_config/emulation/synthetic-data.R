rm(list=ls());gc()

# Load Libraries ####
library(data.table)
library(MicroMoB)
library(ggplot2)
library(progress)

# --------------------------------------------------------------------------------
#   parameters
# --------------------------------------------------------------------------------

# malaria prevalence (from smoothed prevalence data)
pfpr_cluster_mean <- c(0.11315030, 0.07978468, 0.12464991, 0.14256005)

# population in each "patch" (from census data)
N_cluster <- c(182390, 34416, 6455, 2286)

# movement matrtix describing human movement between patches (from travel survey data)
TaR_cluster <- matrix(data = c(
  0.995213840462055, 0.00178539986081067, 0.00272411746835119, 0.000276642208783611,
  0.0134930981974634, 0.984695751678711, 0.00163508718344334, 0.000176062940382268,
  0.0321400695667832, 0.00545417127278543, 0.961687778749627, 0.000717980410804257,
  0.0311841386295664, 0.00685820037404813, 0.00336285202561227, 0.958594808970773
), nrow = 4, ncol = 4, byrow = TRUE)

# movement matrix describing mosquito movement between patches (in general no data exists which can inform this)
psi_cluster <- diag(4)

# mosquito parameters
f <- 0.35 # [0.25, 0.5]
q <- 0.9 # [0.5, 1]
a <- f*q
p = 0.9 # [0.8, 0.95]
g = 1 - p
eip <- 11 # [10-12]
surv <- p^eip

# transmission parameters
b <- 0.55 # [0.45, 0.65]
c <- 0.15 # [0.1, 0.2]

# human infection parameters
r <- 1/200 # [1/250, 1/150]
eta <- 1/32 # [1/40, 1/25]
rho <- 0.06720343 # [0.01, 0.2]

tmax <- 365 * 20


# --------------------------------------------------------------------------------
# initial conditions for model
# --------------------------------------------------------------------------------

pfpr <- pfpr_cluster_mean

odds <- r/(1-rho)*pfpr/(1-(1+rho*r/eta/(1-rho))*pfpr)

# Force of Infection (FOI) in each patch (not on each strata)
h_patch <- as.vector(solve(TaR_cluster) %*% odds)

# EIR in patches
H <- as.vector(N_cluster %*% TaR_cluster)
Z <- (h_patch * H) / (a * b)

I <- pfpr * N_cluster
S <- N_cluster - I - ((rho*I*r)/(eta*(1-rho)))
P <- N_cluster - I - S

# kappa
I_ambient <- as.vector(I %*% TaR_cluster)
kappa <- (I_ambient / H) * c

Y <- Z / surv
M <- (Z*(g + (f*q*p*kappa))) / (f*q*kappa*p*surv) # M from kappa

lambda_eq <- g*M


bt <- function(t,kb,lambdab,phib) {
  kb*exp(-lambdab*((cos(((pi*t) + phib) / (365) ) )^2))
}

kb <- 1
phib <- 50
lambdab <- 5
# plot(bt(t = 1:(365*2),kb = kb,lambdab = lambdab,phib = phib),type='l')

kb_eq <- sapply(X = 1:4, FUN = function(i) {
  optimize(f = function(x, i) {
    (lambda_eq[i] - mean(bt(t = 1:365, kb = x, lambdab = lambdab, phib = phib)))^2
  }, interval = c(1, 1e6), i = i)$minimum
})


lambda <- sapply(X = kb_eq, FUN = function(kb) {
  bt(t = 1:tmax, kb = kb, lambdab = lambdab, phib = phib)
})
lambda <- t(lambda)


# --------------------------------------------------------------------------------
#   simulation
# --------------------------------------------------------------------------------

# initial conditions
n_patch <- nrow(TaR_cluster)
human_state <- matrix(data = 0, nrow = n_patch, ncol = 3, dimnames = list(NULL, c('S', 'I', 'P')))
human_state[, 'S'] <- S
human_state[, 'I'] <- I
human_state[, 'P'] <- P

mod <- make_MicroMoB(tmax = tmax, p = n_patch)

setup_aqua_trace(model = mod, lambda = lambda, stochastic = FALSE)
setup_mosquito_RM(model = mod, stochastic = FALSE, f = f, q = q, eip = eip, p = p, psi = diag(n_patch), nu = 25, M = M, Y = Y, Z = Z)
setup_humans_SIP(model = mod, stochastic = FALSE, theta = TaR_cluster, SIP = human_state, b = b, c = c, r = r, rho = rho, eta = eta)

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
  
  step_aqua(model = mod)
  step_mosquitoes(model = mod)
  step_humans(model = mod)
  
  human_out[day==get_tnow(mod) & state=='S', value := mod$human$SIP[, 'S']]
  human_out[day==get_tnow(mod) & state=='I', value := mod$human$SIP[, 'I']]
  human_out[day==get_tnow(mod) & state=='P', value := mod$human$SIP[, 'P']]
  
  mosy_out[day == get_tnow(mod) & state == 'M', value := mod$mosquito$M]
  mosy_out[day == get_tnow(mod) & state == 'Y', value := mod$mosquito$Y]
  mosy_out[day == get_tnow(mod) & state == 'Z', value := mod$mosquito$Z]
  
  pb$tick()
  mod$global$tnow <- mod$global$tnow + 1L
}

# full output
human_out[, species := "human"]
mosy_out[, species := "mosquito"]
det_out <- rbind(human_out, mosy_out)

# reduced output
human_pfpr <- human_out[, .('pfpr' = .SD[state == 'I', value] / .SD[, sum(value)]), by = .(day, patch)]
mosy_spz <- mosy_out[, .('spz' = .SD[state == 'Z', value] / .SD[state == 'M', value]), by = .(day, patch)]

# diagnostic plots
ggplot(human_pfpr) +
  geom_line(aes(x = day, y = pfpr)) +
  facet_wrap(. ~ patch)

ggplot(mosy_spz) +
  geom_line(aes(x = day, y = spz)) +
  facet_wrap(. ~ patch)


# subsample PfPR
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

sample_beta <- function(mu, var) {
  pars <- estBetaParams(mu, var)
  rbeta(n = length(mu), shape1 = pars[[1]], shape2 = pars[[2]])
}

human_pfpr_samp <- human_pfpr[day >= tmax - (365*2), ]
human_pfpr_samp <- human_pfpr_samp[day %% 20 == 0, ]
human_pfpr_samp[, 'pfpr_samp' := sample_beta(mu = pfpr, var = pfpr * 0.005)]

ggplot(human_pfpr_samp) +
  geom_line(aes(x = day, y = pfpr)) +
  geom_point(aes(x = day, y = pfpr_samp), color = 'red') +
  facet_wrap(. ~ patch, scales = 'free')

# subsample spz
mosy_spz_samp <- mosy_spz[day >= tmax - (365*2), ]
mosy_spz_samp <- mosy_spz_samp[day %% 20 == 0, ]
mosy_spz_samp[, 'spz_samp' := sample_beta(mu = spz, var = spz * 0.005)]

ggplot(mosy_spz_samp) +
  geom_line(aes(x = day, y = spz)) +
  geom_point(aes(x = day, y = spz_samp), color = 'red') +
  facet_wrap(. ~ patch, scales = 'free')

# write data
human_pfpr_csv <- copy(human_pfpr_samp)
human_pfpr_csv[, pfpr := NULL]
human_pfpr_csv[, day := ((0:(.N-1))*20)+1, by = .(patch)]

write.csv2(x = human_pfpr_csv, file = "MicroMoB_config/emulation/synthetic-pfpr.csv")

mosy_spz_csv <- copy(mosy_spz_samp)
mosy_spz_csv[, spz := NULL]
mosy_spz_csv[, day := ((0:(.N-1))*20)+1, by = .(patch)]

write.csv2(x = mosy_spz_csv, file = "MicroMoB_config/emulation/synthetic-spz.csv")


