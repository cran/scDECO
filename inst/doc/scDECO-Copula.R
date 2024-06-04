## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(scDECO)

## -----------------------------------------------------------------------------
n <- 2500

x.use <- rnorm(n)
w.use <- runif(n,-1,1)
marginals.use <- c("ZINB", "ZIGA")

# simulate data
y.use <- scdeco.sim.cop(marginals=marginals.use, x=x.use,
                    eta1.true=c(-2, 0.8), eta2.true=c(-2, 0.8),
                    beta1.true=c(1, 0.5), beta2.true=c(1, 1),
                    alpha1.true=7, alpha2.true=3,
                    tau.true=c(-0.2, .3), w=w.use)

## -----------------------------------------------------------------------------
# fit the model
mcmc.out <- scdeco.cop(y=y.use, x=x.use, marginals=marginals.use, w=w.use,
                       n.mcmc=10, burn=0, thin=1) # n.mcmc=5000, burn=1000, thin=10)

## -----------------------------------------------------------------------------
# extract estimates and confidence intervals
lowerupper <- t(apply(mcmc.out, 2, quantile, c(0.025, 0.5, 0.975)))
estmat <- cbind(lowerupper[,1],
                c(c(-2, 0.8), c(-2, 0.8), c(1, 0.5), c(1, 1), 7, 3, c(-0.2, .3)),
                lowerupper[,c(2,3)])
colnames(estmat) <- c("lower", "trueval", "estval", "upper")
estmat

