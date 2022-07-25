rm(list=ls());gc()

# Load Libraries ####
library(data.table)
library(MicroMoB)


# --------------------------------------------------------------------------------
#   parameters
# --------------------------------------------------------------------------------

# malaria prevalence (from smoothed prevalence data)
pfpr_cluster_mean <- c(0.11315030, 0.07978468, 0.12464991, 0.14256005)
pfpr_cluster_lo <- c(0.0505460, 0.0169550, 0.0259320, 0.0542125)
pfpr_cluster_hi <- c(0.1838235, 0.1746805, 0.2263740, 0.2147102)

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
f <- 0.3 # [0.25, 0.5]
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
lambda <- g*M


# --------------------------------------------------------------------------------
#   simulation
# --------------------------------------------------------------------------------

# initial conditions
n_patch <- nrow(TaR_cluster)
human_state <- matrix(data = 0, nrow = n_patch, ncol = 3, dimnames = list(NULL, c('S', 'I', 'P')))
human_state[, 'S'] <- S
human_state[, 'I'] <- I
human_state[, 'P'] <- P

tmax <- 365 * 3
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
  
  mod$global$tnow <- mod$global$tnow + 1L
}

# full output
human_out[, species := "human"]
mosy_out[, species := "mosquito"]
det_out <- rbind(human_out, mosy_out)

# reduced output
human_pfpr <- human_out[, .('pfpr' = .SD[state == 'I', value] / .SD[, sum(value)]), by = .(day, patch)]
mosy_out <- mosy_out[, .('spz' = .SD[state == 'Z', value] / .SD[state == 'M', value]), by = .(day, patch)]

# diagnostic plots
ggplot(human_pfpr) +
  geom_line(aes(x = day, y = pfpr)) +
  facet_wrap(. ~ patch)

ggplot(mosy_out) +
  geom_line(aes(x = day, y = spz)) +
  facet_wrap(. ~ patch)
