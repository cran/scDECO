## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(scDECO)

## -----------------------------------------------------------------------------
n <- 2500
b.use <- c(-3,0.1)

# simulate the data
simdat <- scdeco.sim.pg(N=n, b0=b.use[1], b1=b.use[2],
                        phi1=4, phi2=4, phi3=1/7,
                        mu1=15, mu2=15, mu3=7,
                        tau0=-2, tau1=0.4)

## -----------------------------------------------------------------------------
# fit the model
mcmc.out <- scdeco.pg(dat=simdat,
                      b0=b.use[1], b1=b.use[2],
                      adapt_iter=1,# 500,
                      update_iter=1, # 500,
                      coda_iter=10, # 5000,
                      coda_thin=1, # 10,
                      coda_burnin=0)# 1000)

## -----------------------------------------------------------------------------
boundsmat <- cbind(mcmc.out$quantiles[,1],
                  c(1/4, 1/4, 7, 15, 15, 7, -2, 0.4), 
                  mcmc.out$quantiles[,c(3,5)])

colnames(boundsmat) <- c("lower", "true", "est", "upper")

boundsmat

