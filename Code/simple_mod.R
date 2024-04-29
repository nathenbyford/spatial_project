## Nathen Byford
## Simple poisson regression ignoring spatial component

library(tidyverse) # For reproducing plots seen in the paper.
library(nimble) # For MCMC computation using NIMBLE.
library(coda) # For manipulation of MCMC results.
library(mgcv)
library(ngspatial)
library(sp) # For plotting the micro-regions.
library(spdep) # For computing the neighbourhood and adjancency objects.
library(maps) # For adding a map to the plots of Brazil.
library(viridis)
library(mapproj)

set.seed(794637) # Set the seed (found in the Master script).

load("Code/covid.RData") # Load in the data. STILL NEED TO CHANGE THIS
load("spatial.Rdata")

# Define the model code
code <- nimbleCode({
  # Likelihood: Poisson model
  for (i in 1:n) {
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha + beta1 * c1[i] + beta2 * c2[i]
  }
  
  # Priors
  alpha ~ dnorm(0, sd = 10)  # Normal prior for intercept
  beta1 ~ dnorm(0, sd = 10)  # Normal prior for slope
  beta2 ~ dnorm(0, sd = 10)  # Normal prior for slope
})

constants <- list(n = 49, c1 = dat$Excessive_Drinking,
                  c2 = States$Unemployment)
data <- list(y = dat$positive)

inits <- list(alpha = 0, beta1 = 0, beta2 = 0)

model <- nimbleModel(code, data = data, constants = constants, inits = inits)
mod_conf <- configureMCMC(model, enableWAIC = TRUE)

mcmc <- buildMCMC(mod_conf)

compiled_mod <- compileNimble(model)
compiled_mcmc <- compileNimble(mcmc, project = model)

samples <- runMCMC(compiled_mcmc, nburnin = 2000, niter = 5000,
                  nchains = 3, WAIC = TRUE)

samples
