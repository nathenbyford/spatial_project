library(R2OpenBUGS)
library(coda)

# Load data

load("covid.RData")
load("spatial.Rdata")

# Set up necessary data for the ICAR prior.
adjacency=unlist(neighborhood)
n_adj=card(neighborhood)
l_adj=length(adjacency)
weights=rep(1,l_adj)

Covid_N=length(dat$positive) # Number of observations.
n_regions=length(n_adj) # Number of regions.

total_obs=c(sum(TBdata$TB[1:n_regions]),sum(TBdata$TB[(1:n_regions)+n_regions]),sum(TBdata$TB[(1:n_regions)+2*n_regions]))

# Set up index for spatial parameters rho and delta.
region_index=numeric(n_regions)
for(i in 1:n_regions){
  region_index[i]=which(States$NAME_1==dat$NAME_1[i])
}
region_index=rep(region_index,3)


# Create polynomials.
poly_unem=poly(States$Unemployment,2)
poly_den=poly(dat$Popdensity,2)


model_data <- list(
  z = States$Covid.2020,
  n = Covid_N, pop = dat$Pop, n_adj = n_adj, adj = adjacency,
  unem = poly_unem, den = poly_den, R = n_regions,
  index = region_index, l_adj = length(adjacency)
)

bugs_code <- "model {
  for (i in 1:n) {
    mu_pi[i] <- b[1] + gamma[i]
    pi[i] <- ilogit(mu_pi[i])

    lambda[i] <- exp(log(pop[i]) + a[1] + unem[i,1]*a[2] + unem[i,2]*a[3] +
                       den[i,1]*a[4] + den[i,2]*a[5] + phi[index[i]] + theta[index[i]])

    z[i] ~ dpois(pi[i]*lambda[i])
    gamma[i] ~ dnorm(0, 1/epsilon^2)
  }

  for (j in 1:R) {
    theta[j] ~ dnorm(0, 1/sigma^2)
  }

  for (j in 1:R) {
    phi[j] ~ car.normal(adj[1:l_adj], num[1:R], tau)
  }

  a[1] ~ dnorm(-8, 1)
  for (i in 2:6) {
    a[i] ~ dnorm(0, 10)
  }
  
  b[1] ~ dnorm(2, 0.6)
  
  sigma ~ dgamma(0.5, 0.0005)  # precision parameter for theta
  epsilon ~ dgamma(0.5, 0.0005) # precision parameter for gamma
  nu ~ dgamma(0.5, 0.0005)      # precision parameter for phi
  tau <- 1/nu^2
}
"

cat(bugs_code, file = "Covid_bugs.txt")

Covid_inits1=list(sigma=0.25,nu=1,epsilon=0.25,a=c(-7,rep(-0.1,5)),
                  gamma=rnorm(Covid_N,0,2.5),phi=rnorm(n_regions,0,10),theta=rnorm(n_regions,0,5))
Covid_inits2=list(sigma=0.5,nu=0.75,epsilon=0.5,a=c(-7,rep(0.1,5)),
                  gamma=rnorm(Covid_N,0,5),phi=rnorm(n_regions,0,7.5),theta=rnorm(n_regions,0,10))
Covid_inits3=list(sigma=0.75,nu=0.5,epsilon=0.25,a=c(-9,rep(-0.1,5)),
                  gamma=rnorm(Covid_N,0,2.5),phi=rnorm(n_regions,0,5),theta=rnorm(n_regions,0,15))
Covid_inits4=list(sigma=1,nu=0.25,epsilon=0.5,a=c(-9,rep(0.1,5)),
                  gamma=rnorm(Covid_N,0,5),phi=rnorm(n_regions,0,2.5),theta=rnorm(n_regions,0,20))
Covid_inits=list(Covid_inits1, Covid_inits2, Covid_inits3, Covid_inits4)

parameters <- c("gamma", "theta", "phi", "a", "b", "sigma", "epsilon", "nu")

model <- bugs(data = model_data,
              inits = Covid_inits,
              parameters.to.save = parameters,
              model.file = "Covid_bugs.txt",
              n.chains = 4, n.iter = 20000,
              n.burnin = 10000, debug = TRUE)
