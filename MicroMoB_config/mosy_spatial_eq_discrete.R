rm(list=ls()); gc()
library(deSolve)
library(expm)
library(MASS)
library(data.table)
library(ggplot2)

p <- 1-(1/20)
g <- 1 - p

n_patch <- 4

# movement matrix
K <- matrix(0, n_patch, n_patch)
K[upper.tri(K)] <- rexp(sum(1:(n_patch-1)))
K[lower.tri(K)] <- rexp(sum(1:(n_patch-1)))
diag(K) <- 10
K <- K/rowSums(K)

M <- rpois(n = n_patch, lambda = 10)

M %*% K
sum(M)
sum(M %*% K)


