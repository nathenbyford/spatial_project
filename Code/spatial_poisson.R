## Nathen Byford
## Spatial poisson regression model for project

library(tidyverse) # For reproducing plots seen in the paper.
library(nimble) # For MCMC computation using NIMBLE.
library(coda) # For manipulation of MCMC results.
library(mgcv)
library(ngspatial)
library(spdep) # For computing the neighbourhood and adjancency objects.
library(maps) # For adding a map to the plots of Brazil.
library(viridis)
library(mapproj)

load("Code/covid.RData")
load("data/spatial.Rdata")

set.seed(794637)

n <- 49 # number of observations
adjacency <- unlist(neighborhood)
adj_matrix <- nb2mat(neighborhood, style = "B")
weight_matrix <- adj_matrix
n_adj <- card(neighborhood)
l_adj <- length(adjacency)
weights <- rep(1,l_adj)

n_regions <- length(n_adj)


# Model code for spatial poisson regression
code <- nimbleCode({
  # Likelihood: Poisson model
  for (i in 1:n) {
    counts[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha + beta * covariate[i] + spatial_effect[i]
  }
  
  # Priors
  alpha ~ dnorm(0, sd = 1000)  # Intercept prior
  beta ~ dnorm(0, sd = 1000)   # Covariate effect prior
  
  # Spatial effect prior (Conditional autoregressive prior)
  for (i in 1:n) {
    spatial_effect[i] ~ dcar_normal(
      adj = adj_matrix[i, 1:neighbor_count[i]],  # Neighbors of state i
      weights = weights_matrix[i, 1:neighbor_count[i]],  # Weights for neighbors
      num = neighbor_count[i],  # Number of neighbors
      tau = tau
    )
  }
  tau ~ dgamma(0.5, 0.0005)  # Precision parameter for spatial effect
  
})

adj_list <- apply(adj_matrix, 2, \(x) x[x != 0])
weight_list <- apply(weight_matrix, 2, \(x) as.numeric(x[x != 0]))
neighbor_count <- rowSums(adj_matrix != 0.0)

# Specify constants
constants <- list(n = n,
                  neighbor_count = n_adj,
                  adj_matrix = adj_matrix,
                  weights_matrix = weight_matrix)
data <- list(counts = dat$positive, covariate = dat$Popdensity)

# Model
model <- nimbleModel(code = code,
                     constants = constants,
                     data = data)
c_model <- compileNimble(model)



## Trying in stan

library(rstan)

stan_code <- "
// function of ICAR following Morris~et~al.
functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1,
  int[] node2) {
    return−0.5 * dot_self(phi[node1] − phi[node2]) +
    normal_lpdf(sum(phi) | 0,0.001 * N);
  }
}

data {
    int<lower=1> n;  // Number of regions
    int counts[n];  // Observed counts in each region
    vector[n] covariate;  // Covariate data for each region
    int<lower=0, upper=n> neighbor_count[n];  // Number of neighbors for each region
    int adj[n, max(neighbor_count)];  // Adjacency list for each region
    vector[n] weights[n];  // Weights for each region's neighbors
}

parameters {
    real alpha;  // Intercept
    real beta;  // Coefficient for covariate
    vector[n] spatial_effect;  // Spatial effects
    real<lower=0> tau;  // Precision parameter for CAR prior
}

model {
    // Likelihood: Poisson model
    for (i in 1:n) {
        counts[i] ~ poisson_log(alpha + beta * covariate[i] + spatial_effect[i]);
    }

    // Priors
    alpha ~ normal(0, 1000);  // Prior for intercept
    beta ~ normal(0, 1000);  // Prior for coefficient

    // Spatial effect prior using CAR distribution
    spatial_effect ~ icar_normal(adj, weights, tau, neighbor_count);

    // Prior for precision parameter
    tau ~ gamma(0.5, 0.5);
}
"