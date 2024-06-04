#' Simulating from copula model
#'
#' @param marginals provide vector of length 2 of which marginals to use
#' @param x covariate matrix
#' @param eta1.true zero-inflation parameters for marginal 1
#' @param eta2.true zero-inflation parameters for marginal 2
#' @param beta1.true mean coefficients for marginal 1
#' @param beta2.true mean coefficients for marginal 2
#' @param alpha1.true second parameter coefficients for marginal 1
#' @param alpha2.true second parameter coefficients for marginal 2
#' @param tau.true coefficients for correlation
#' @param w (optional) covariate matrix for zero-inflation portion
#'
#' @return matrix with values simulated from copula model
#' @export
#'
#' @examples
#' n <- 2500
#' x.use = rnorm(n)
#' w.use = runif(n,-1,1)
#' eta1.use = c(-2.2, 0.7)
#' eta2.use = c(-2, 0.8)
#' beta1.use = c(1,0.5)
#' beta2.use = c(1,1)
#' alpha1.use = 7
#' alpha2.use = 3
#' tau.use = c(-0.2, .3)
#'
#' marginals.use <- c("ZINB", "ZIGA")
#'
#' y.use <- scdeco.sim.cop(marginals=marginals.use, x=x.use,
#'                     eta1.true=eta1.use, eta2.true=eta2.use,
#'                     beta1.true=beta1.use, beta2.true=beta2.use,
#'                     alpha1.true=alpha1.use, alpha2.true=alpha2.use,
#'                     tau.true=tau.use, w=w.use)
#'
#' y.use[1:10,]
#'
scdeco.sim.cop <- function(marginals, x, eta1.true, eta2.true, beta1.true, beta2.true, alpha1.true, alpha2.true, tau.true, w=NULL){

  x <- cbind(1, x)
  n <- nrow(x)

  if ((length(w)==0) && (is.null(w))){
    w <- runif(n)
  }

  if (NROW(w) != n){
    stop("w must have same NROW as y")
  }

  w <- cbind(1, w)

  qmarg1 <- qmarg_finder(marginals[1])
  qmarg2 <- qmarg_finder(marginals[2])


  meanlink1 <- meanlink_finder(marginals[1])
  meanlink2 <- meanlink_finder(marginals[2])

  mu1 = c(meanlink1(x%*%beta1.true))
  mu2 = c(meanlink2(x%*%beta2.true))
  p1 = c(1/(exp(-w%*%eta1.true)+1))
  p2 = c(1/(exp(-w%*%eta2.true)+1))

  if (!ZImarg(marginals[1])){
    p1 <- p1*0
  }

  if (!ZImarg(marginals[2])){
    p2 <- p2*0
  }

  rho = c(exp(x%*%tau.true)-1)/
    c(exp(x%*%tau.true)+1)


  y = matrix(0,n,2)
  for(i in 1:n){
    R = matrix(c(1,rho[i],rho[i],1),2,2)
    z = c(rnorm(2)%*%chol(R))

    u = pnorm(z)

    y[i,1] = ifelse(u[1] < p1[i], 0,
                    qmarg1((u[1]-p1[i])/(1-p1[i]), mu1[i], alpha1.true))

    y[i,2] = ifelse(u[2] < p2[i], 0,
                    qmarg2((u[2]-p2[i])/(1-p2[i]), mu2[i], alpha2.true))


  }

  if (marginals[1] %in% c("Beta", "ZIBeta")){
    y[,1] <- beta_fixer(y[,1])
  }

  if (marginals[2] %in% c("Beta", "ZIBeta")){
    y[,2] <- beta_fixer(y[,2])
  }


  return(y)

}





beta_fixer <- function(y){
  y[y>0.9999] <- y[y>0.9999] - 0.001
  return(y)
}
