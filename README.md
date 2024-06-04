<p align="center">
  <img src="https://raw.githubusercontent.com/abussing/abussing/abb5400cf82e44abc85d612a10d9ec8a2ed65e09/images/scdeco_logo.svg" alt="scDECO logo" width="400">
</p>

<p align="center">
  <align="center"><ins>s</ins>ingle-<ins>c</ins>ell <ins>d</ins>iff<ins>e</ins>rential <ins>co</ins>-expression
</p>

<div align="center">
  <a href="https://doi.org/10.1111/biom.13701">
    <img src="https://img.shields.io/badge/DOI-doi.org%2F10.1111%2Fbiom.13701-garnet?color=%23840028" alt="DOI">
  </a>
  <a href="https://doi.org/10.1111/biom.13457">
    <img src="https://img.shields.io/badge/DOI-doi.org%2F10.1111%2Fbiom.13457-garnet?color=%23840028" alt="DOI">
  </a>
</div>





## Description

scDECO is an R package for estimating differential co-expression in single-cell RNA-seq data. 

Differential co-expression refers to correlation between two random variables $X_1, X_2$ changing with respect to the value of some covariate(s). 

<p align="center">
<img src="https://raw.githubusercontent.com/abussing/abussing/abb5400cf82e44abc85d612a10d9ec8a2ed65e09/images/dynamic_corr_plot.svg" alt="differential coexpression" width="700">
</p>

The package contains implementations for two different Bayesian models:
1. `scDECO.cop`: fits a bivariate Gaussian copula model with flexible, covariate-dependent, optionally zero-inflated marginals
2. `scDECO.pg`: fits a zero-inflated bivariate Poisson-Gamma model with correlation imparted through a latent bivariate normal variable

## Tutorial

Quick video tutorial on the package:

[![youtube image](https://img.youtube.com/vi/tsmDzepHIVg/0.jpg)](https://www.youtube.com/watch?v=tsmDzepHIVg)






## Installation

```{r}
devtools::install_github("YenYiHo-Lab/scDECO")
library(scDECO)
```

## Usage

### `scdeco.cop`

```{r}
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
```
Parameters:
* `marginals`: The two marginals. Options are NB, ZINB, GA, ZIGA, Beta, ZIBEta
* `x`: The vector (or matrix) containing the covariate values to be regressed for mean and rho parameters.
* `eta1.true`: The coefficients of the 1st marginal's zero-inflation parameter. 
* `eta2.true`: The coefficients of the 2nd marginal's zero-inflation parameter. 
* `beta1.true`: The coefficients of the 1st marginal's mean parameter. 
* `beta2.true`: The coefficients of the 2nd marginal's mean parameter. 
* `alpha1.true`: The coefficient of the 1st marginal's second parameter. 
* `alpha2.true`: The coefficient of the 2nd marginal's second parameter.
* `tau.true`: The coefficients of the correlation parameter. 
* `w`: A vector (or matrix) containing the covariate values to be regressed for zero-inflation parameters.

This will simulate a 2-column matrix of `NROW(x)` rows of observations from the scdeco.cop model.


```{r}
# fit the model
mcmc.out <- scdeco.cop(y=y.use, x=x.use, marginals=marginals.use, w=w.use,
                     n.mcmc=5000, burn=1000, thin=10)
```
Parameters:
* `y`: 2-column matrix with the dependent variable observations.
* `n.mcmc`: The number of MCMC iterations to run.
* `burn`: The number of MCMC iterations to burn from the beginning of the chain.
* `thin`: The number of MCMC iterations to thin.

This will return a matrix where the columns correspond to the different parameters of the model and the rows correspond to MCMC samples where the burn and thin has already been incorporated. 

One can obtain estimates and confidence intervals for each parameter by looking at quantiles of these MCMC samples.

```{r}
# extract estimates and confidence intervals
lowerupper <- t(apply(mcmc.out, 2, quantile, c(0.025, 0.5, 0.975)))
estmat <- cbind(lowerupper[,1],
                c(c(-2, 0.8), c(-2, 0.8), c(1, 0.5), c(1, 1), 7, 3, c(-0.2, .3)),
                lowerupper[,c(2,3)])
colnames(estmat) <- c("lower", "trueval", "estval", "upper")
estmat
```

### `scdeco.pg`

```{r}
n <- 2500
b.use <- c(-3,0.1)

# simulate the data
simdat <- scdeco.sim.pg(N=n, b0=b.use[1], b1=b.use[2],
                        phi1=4, phi2=4, phi3=1/7,
                        mu1=15, mu2=15, mu3=7,
                        tau0=-2, tau1=0.4)
```
Parameters:
* `N`: Sample size for the simulated data.
* `b0`: The intercept coefficient of the zero-inflation parameter. 
* `b1`: The slope coefficient of the zero-inflation parameter. 
* `phi1`: The over-dispersion parameter of the 1st ZINB marginal.
* `phi2`: The over-dispersion parameter of the 2nd ZINB marginal.
* `phi3`: The over-dispersion parameter of the ZINB covariate vector. 
* `mu1`: The mean parameter of the 1st ZINB marginal.
* `mu2`: The mean parameter of the 2nd ZINB marginal.
* `mu3`: The mean parameter of the ZINB covariate vector.
* `tau0`: The intercept coefficient of the correlation parameter. 
* `tau1`: The slope coefficient of the correlation parameter.

This will simulate a 3-column matrix of $N$ rows, where the first two columns are observations and the third column is the ZINB covariate which will be used in regressing the correlation parameter of the scdeco.pg model. 

```{r}
# fit the model
mcmc.out <- scdeco.pg(dat=simdat,
                      b0=b.use[1], b1=b.use[2],
                      adapt_iter=500,
                      update_iter=500,
                      coda_iter=5000,
                      coda_thin=10,
                      coda_burnin=1000)
```
Parameters:
* `dat`: The 3-column matrix where the first two columns are observations and the third column is the ZINB covariate. An additional covariate can be added as a 4th column if desired.
* `adapt_iter`: The number of adaptive iterations to run.
* `update_iter`: The number of update iterations to run.
* `coda_iter`: The number of MCMC iterations to run after the adapt and update.
* `coda_thin`: The number of MCMC iterations to burn from the coda_iter iterations.
* `coda_burnin`: The number of MCMC iterations to thin from the coda_burnin iterations.
  
This will return a matrix where the columns correspond to the different parameters of the model and the rows correspond to MCMC samples where the adapt, update, burn, and thin has already been incorporated. 

One can obtain estimates and confidence intervals for each parameter by looking at quantiles of these MCMC samples using the same code as in the scdeco.cop example.



## Model Details

For details on the exact models, please see the vignettes of the scDECO R package.

## Please Cite

Zichen Ma, Shannon W. Davis, Yen-Yi Ho, Flexible Copula Model for Integrating Correlated Multi-Omics Data from Single-Cell Experiments, Biometrics, Volume 79, Issue 2, June 2023, Pages 1559–1572, https://doi.org/10.1111/biom.13701

Zhen Yang, Yen-Yi Ho, Modeling Dynamic Correlation in Zero-Inflated Bivariate Count Data with Applications to Single-Cell RNA Sequencing Data, Biometrics, Volume 78, Issue 2, June 2022, Pages 766–776, https://doi.org/10.1111/biom.13457

## Contact

If you have any questions or found any issues, please contact: abussing@email.sc.edu.





