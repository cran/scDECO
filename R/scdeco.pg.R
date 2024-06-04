#' ZENCO Poisson Gamma dynamic correlation fitting function
#'
#' @param dat matrix containing expression values as first two columns and covariate as third column
#' @param b0 intercept of zinf parameter
#' @param b1 slope of zinf parameter
#' @param adapt_iter number of adaptation iterations in the jags.model function
#' @param update_iter update iterations in the update function
#' @param coda_iter number of iterations for the coda.sample function
#' @param coda_thin how much to thin the resulting MCMC output
#' @param coda_burnin how many iterations to burn before beginning coda sample collection
#'
#' @return MCMC samples that have been adapted, burned, and thinned
#' @export
#'
#' @examples
#'
#' phi1_use <- 4
#' phi2_use <- 4
#' phi3_use <- 1/7
#' mu1_use <- 15
#' mu2_use <- 15
#' mu3_use <- 7
#' b0_use <- -3
#' b1_use <- 0.1
#' tau0_use <- -2
#' tau1_use <- 0.4
#'
#' simdat <- scdeco.sim.pg(N=1000, b0=b0_use, b1=b1_use,
#'                         phi1=phi1_use, phi2=phi2_use, phi3=phi3_use,
#'                         mu1=mu1_use, mu2=mu2_use, mu3=mu3_use,
#'                         tau0=tau0_use, tau1=tau1_use)
#'
#' zenco_out <- scdeco.pg(dat=simdat,
#'                        b0=b0_use, b1=b1_use,
#'                        adapt_iter=1, # 500,
#'                        update_iter=1, # 500,
#'                        coda_iter=5, # 5000,
#'                        coda_thin=1, # 10,
#'                        coda_burnin=0) # 1000
#'
#' boundsmat <- cbind(zenco_out$quantiles[,1],
#'                    c(1/phi1_use, 1/phi2_use, 1/phi3_use,
#'                    mu1_use, mu2_use, mu3_use,
#'                    tau0_use, tau1_use),
#'                    zenco_out$quantiles[,c(3,5)])
#'
#' colnames(boundsmat) <- c("lower", "true", "est", "upper")
#'
#' boundsmat
#'
scdeco.pg <- function(dat, b0, b1, adapt_iter=100, update_iter=100, coda_iter=1000, coda_thin=5, coda_burnin=100){
  if (ncol(dat) == 4){
    scdeco.pg.xc(dat, b0, b1, adapt_iter, update_iter, coda_iter, coda_thin, coda_burnin)
  } else if (ncol(dat) == 3){
    scdeco.pg.noxc(dat, b0, b1, adapt_iter, update_iter, coda_iter, coda_thin, coda_burnin)
  } else {
    stop("dat needs 3 or 4 columns")
  }
}
