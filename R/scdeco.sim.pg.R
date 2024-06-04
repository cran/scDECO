#' Simulating from ZENCO Model
#'
#' @param N size of sample to be generated
#' @param b0 intercept of zinf parameter
#' @param b1 slope of zinf parameter
#' @param phi1 over-dispersion parameter of first marginal
#' @param phi2 over-dispersion parameter of second marginal
#' @param phi3 over-dispersion parameter of covariate vector
#' @param mu1 mean parameter of first marginal
#' @param mu2 mean parameter of second marginal
#' @param mu3 mean parameter of covariate vector
#' @param tau0 intercept of correlation
#' @param tau1 slope of of correlation
#' @param tau2 (optional) correlation coefficient on optional xc covariate vector
#' @param tau3 (optional) correlation coefficient on interaction between x3 and xc
#' @param xc (optional) secondary covariate to be regressed
#'
#' @return a matrix with expressions as first two columns and covariates as remaining columns
#' @export
#' @import MASS
#'
#' @examples
#' phi1_use <- 4
#' phi2_use <- 4
#' phi3_use <- 1/6
#' mu1_use <- 15
#' mu2_use <- 15
#' mu3_use <- 7
#' b0_use <- 0.6882
#' b1_use <- -0.2995
#' tau0_use <- 0.07
#' tau1_use <- 0.05
#'
#' simdat <- scdeco.sim.pg(N=1000, b0=b0_use, b1=b1_use,
#'                         phi1=phi1_use, phi2=phi2_use, phi3=phi3_use,
#'                         mu1=mu1_use, mu2=mu2_use, mu3=mu3_use,
#'                         tau0=tau0_use, tau1=tau1_use)
#' simdat[1:10,]
scdeco.sim.pg <- function(N, b0, b1, phi1, phi2, mu1, mu2, tau0, tau1, mu3, phi3, tau2=NULL, tau3=NULL, xc=NULL){

  if ((length(xc) > 0) & (length(xc) != N)){
    stop("xc must be same length as N")
  }

  phi <- c(phi1,phi2)
  mu <- c(mu1,mu2)
  lambda0 <- 1e-10 # mean of poisson
  # dropout rate
  p <- c(exp(b0+b1*mu[1])/(1+exp(b0+b1*mu[1])), exp(b0+b1*mu[2])/(1+exp(b0+b1*mu[2])))
  cp <- mu
  # x3
  x3 <- rnbinom(N, size=1/phi3, mu=mu3)
  p3 <- rbinom(N, size=1, prob=exp(b0+b1*mu3)/(1+exp(b0+b1*mu3)))
  z3 <- rpois(N, lambda0) # poisson
  x3 <- (1-p3)*x3+p3*z3

  if (length(xc) == 0){
    # if x3 dropout, x3 is mu3,
    rho.mu3 <- (exp(tau0+tau1*mu3)-1)/(exp(tau0+tau1*mu3)+1)
    rho <- (exp(tau0+tau1*x3)-1)/(exp(tau0+tau1*x3)+1)
  } else {
    mu.c <- mean(xc)
    # if x3 dropout, x3 is mu3, xc is mean of x3
    rho.mu3 <- (exp(tau0+tau1*mu3 + tau2*mu.c + tau3*mu.c*mu3)-1)/(exp(tau0+tau1*mu3 + tau2*mu.c + tau3*mu.c*mu3)+1)
    rho <- (exp(tau0+tau1*x3 + tau2*xc + tau3*xc*x3)-1)/(exp(tau0+tau1*x3 + tau2*xc + tau3*xc*x3)+1)
  }
  rho <- (1-p3)*rho+p3*rho.mu3
  # quadratic
  cp[1] <- 1/phi[1]
  cp[2] <- 1/phi[2]
  z <- matrix(rep(0,2*N), ncol=2)
  y <- matrix(rep(0,2*N), ncol=2)
  for(i in 1:N){
    Sigma <- matrix(c(1,rho[i],rho[i],1),2,2)
    z[i, ] <- mvrnorm(n = 1, rep(0, 2), Sigma)
    for(j in 1:2){
      y[i,j] <- rpois(1,qgamma(pnorm(z[i,j]),cp[j],cp[j])*mu[j])
    }
  }

  # Add "Zero"
  z1 <- rpois(N, lambda0) # poisson
  z2 <- rpois(N, lambda0)
  p1 <- rbinom(N, size=1, prob=p[1])
  p2 <- rbinom(N, size=1, prob=p[2])
  y1 <- (1-p1)*y[,1]+p1*z1
  y2 <- (1-p2)*y[,2]+p2*z2
  y <- cbind(y1, y2)
  return(cbind(y,x3,xc))
}
