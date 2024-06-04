#' Copula dynamic correlation fitting function
#'
#' @param y 2-column matrix of observations
#' @param x covariates
#' @param marginals length-2 vector with strings of the two marginals
#' @param w (optional)
#' @param n.mcmc number of mcmc iterations to run
#' @param burn how many of the mcmc iterations to burn
#' @param thin how much to thin the mcmc iterations
#'
#' @import msm
#'
#' @importFrom stats cor cov dbeta dgamma dnorm pbeta pgamma pnbinom pnorm qbeta qgamma qnbinom qnorm rbinom rnbinom rnorm rpois runif sd update var
#'
#'
#' @return matrix with mcmc samples as rows and columns corresponding to the different parameters
#' @export
#'
#' @examples
#' n <- 1000
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
#' mcmc.out <- scdeco.cop(y=y.use, x=x.use, marginals=marginals.use, w=w.use,
#'                       n.mcmc=10, burn=0, thin=1) # n.mcmc=1000, burn=100, thin=5)
#'
#' lowerupper <- t(apply(mcmc.out, 2, quantile, c(0.025, 0.5, 0.975)))
#' estmat <- cbind(lowerupper[,1],
#'                 c(eta1.use, eta2.use, beta1.use, beta2.use, alpha1.use, alpha2.use, tau.use),
#'                 lowerupper[,c(2,3)])
#' colnames(estmat) <- c("lower", "trueval", "estval", "upper")
#' estmat
#'
scdeco.cop <- function(y,x,marginals,w=NULL, n.mcmc=5000, burn=1000, thin=10){

  n = nrow(y)

  if ((length(w)==0) && (is.null(w))){
    w <- runif(n)
  }

  if (NROW(w) != n){
    stop("w must have same NROW as y")
  }

  w <- cbind(1, w)

  if (all(w[,2] == 1)){
    if (ncol(w)==2){
      stop("w provided had all ones")
    } else {
      w <- w[,-1]
    }
  }


  if (NROW(x) != n){
    stop("x must have same NROW as y")
  }

  x <- cbind(1, x)

  if (all(x[,2] == 1)){
    if (ncol(x)==2){
      stop("x provided had all ones")
    } else {
      x <- x[,-1]
    }
  }


  d = ncol(x)
  d.w <- ncol(w)


  alpha1.update <- alphafun_finder(marginals[1])
  alpha2.update <- alphafun_finder(marginals[2])

  z1.update <- zupdate_finder(marginals[1])
  z2.update <- zupdate_finder(marginals[2])


  beta1 = beta2 = tau = matrix(0,n.mcmc,d)
  alpha1 = alpha2 = rep(1,n.mcmc)
  eta1 = matrix(0,n.mcmc,d.w)
  p1 = c(exp(w%*%(eta1[1,]))/(1+exp(w%*%(eta1[1,]))))
  eta2 = matrix(0,n.mcmc,d.w)
  p2 = c(exp(w%*%(eta2[1,]))/(1+exp(w%*%(eta2[1,]))))

  mu1 = mu2 = c(exp(x%*%(beta1[1,])))
  rho = c((exp(x%*%(tau[1,]))-1)/(1+exp(x%*%(tau[1,]))))
  z = matrix(rnorm(n*2),n,2)%*%chol(cor(y))

  if (!ZImarg(marginals[1])){
    p1 <- p1*0
  }

  if (!ZImarg(marginals[2])){
    p2 <- p2*0
  }


  # Loop n.mcmc times:
  for(i in 2:n.mcmc){
    beta1.mcmc = beta.update(marginals[1],beta1[1:(i-1),],y[,1],x,z[,2],mu1,
                             alpha1[i-1],p1,rho)
    beta1[i,] = beta1.mcmc$beta
    mu1 = beta1.mcmc$mu

    alpha1[i] = alpha1.update(alpha1[1:(i-1)],y[,1],z[,2],mu1,p1,rho)

    eta1.mcmc = eta.update(marginals[1],eta1[1:(i-1),],y[,1],w,z[,2],mu1,
                           alpha1[i],p1,rho)

    eta1[i,] = eta1.mcmc$eta
    p1 = eta1.mcmc$p

    z[,1] = z1.update(y[,1],z[,2],mu1,alpha1[i],p1,rho)

    beta2.mcmc = beta.update(marginals[2],beta2[1:(i-1),],y[,2],x,z[,1],mu2,
                             alpha2[i-1],p2,rho)
    beta2[i,] = beta2.mcmc$beta
    mu2 = beta2.mcmc$mu

    alpha2[i] = alpha2.update(alpha2[1:(i-1)],y[,2],z[,1],mu2,p2,rho)

    eta2.mcmc = eta.update(marginals[2],eta2[1:(i-1),],y[,2],w,z[,1],mu2,
                           alpha2[i],p2,rho)
    eta2[i,] = eta2.mcmc$eta
    p2 = eta2.mcmc$p

    z[,2] = z2.update(y[,2],z[,1],mu2,alpha2[i],p2,rho)

    tau.mcmc = tau.update(tau[1:(i-1),],z,x,rho)
    tau[i,] = tau.mcmc$tau
    rho = tau.mcmc$rho
  }
  # samp = seq(burn+1,n.mcmc,by=thin)
  mcmc.out = cbind(eta1,eta2,beta1,beta2,alpha1,alpha2,tau)#[samp,]


  colnames(mcmc.out) = c(paste0("eta1",0:(ncol(eta1)-1)), paste0("eta2",0:(ncol(eta2)-1)),
                         paste0("beta1",0:(ncol(beta1)-1)), paste0("beta2",0:(ncol(beta2)-1)),
                         "alpha1","alpha2",
                         paste0("tau",0:(ncol(tau)-1)))
  return(mcmc.out)}












pgam = function(y,mu,lambda,p){(1-p)*pgamma(y,mu*lambda,lambda) + p}



pbet = function(y,mu,phi,p){(1-p)*pbeta(y,mu*phi,(1-mu)*phi) + p}



mvdens.z = function(z,rho){
  llik = sum(-.5*(z[,1]^2 - 2*rho*z[,1]*z[,2] + z[,2]^2)/(1-rho^2) -
               .5*log(1-rho^2),na.rm=T)
  return(llik)}



ZImarg <- function(marginal){

  marginal %in% c("ZINB", "ZIGA", "ZIBeta")

}



qmarg_finder <- function(marginal){

  if (marginal %in% c("NB","ZINB")){
    qmarg <- function(x, mu, alpha){qnbinom(p=x, mu=mu, size=alpha)}
  } else if (marginal %in% c("GA","ZIGA")){
    qmarg <- function(x, mu, alpha){qgamma(x,mu*alpha,alpha)}
  } else if (marginal %in% c("Beta","ZIBeta")){
    qmarg <- function(x,mu,alpha){qbeta(x,mu*alpha,(1-mu)*alpha)}
  } else {
    stop("marginal is not one of NB, ZINB, GA, ZIGA, Beta, ZIBeta")
  }

  return(qmarg)

}



llike_finder <- function(marginal){

  if (marginal %in% c("NB","ZINB")){
    llike <- llike.NB
  } else if (marginal %in% c("GA","ZIGA")){
    llike <- llike.GA
  } else if (marginal %in% c("Beta","ZIBeta")){
    llike <- llike.Beta
  } else {
    stop("marginal is not one of NB, ZINB, GA, ZIGA, Beta, ZIBeta")
  }

  return(llike)

}



zupdate_finder <- function(marginal){

  if (marginal %in% c("NB","ZINB")){
    llike <- zNB.update
  } else if (marginal %in% c("GA","ZIGA")){
    llike <- zGA.update
  } else if (marginal %in% c("Beta","ZIBeta")){
    llike <- zBeta.update
  } else {
    stop("marginal is not one of NB, ZINB, GA, ZIGA, Beta, ZIBeta")
  }

  return(llike)

}



alphafun_finder <- function(marginal){

  if (marginal %in% c("NB","ZINB")){
    llike <- alphaNB.update
  } else if (marginal %in% c("GA","ZIGA")){
    llike <- alphaGA.update
  } else if (marginal %in% c("Beta","ZIBeta")){
    llike <- alphaBeta.update
  } else {
    stop("marginal is not one of NB, ZINB, GA, ZIGA, Beta, ZIBeta")
  }

  return(llike)

}




meanlink_finder <- function(marginal){

  if (marginal %in% c("Beta","ZIBeta")){
    meanlink <- function(x){exp(x)/(1+exp(x))}
  } else {
    meanlink <- function(x){exp(x)}
  }

  return(meanlink)

}







llike.NB = function(y,z,mu,alpha,p,rho){
  id1 = which(y!=0)
  id0 = which(y==0)

  uu = pnbinom(y[id1],mu=mu[id1],size=alpha)*(1-p[id1]) + p[id1]
  ul = pnbinom(y[id1]-1,mu=mu[id1],size=alpha)*(1-p[id1]) + p[id1]

  lik = pnorm(qnorm(uu),mean=rho[id1]*z[id1],sd=sqrt(1-rho[id1]^2)) -
    pnorm(qnorm(ul),mean=rho[id1]*z[id1],sd=sqrt(1-rho[id1]^2))
  lik = ifelse(lik==0,.00001,lik)
  llik = sum(log(lik),na.rm=T)

  if(length(id0) > 0){
    u = pnbinom(0,mu=mu[id0],size=alpha)*(1-p[id0]) + p[id0]
    lik = pnorm(qnorm(u),mean=rho[id0]*z[id0],sd=sqrt(1-rho[id0]^2))
    lik = ifelse(lik==0,lik+.00001,lik)
    llik = llik + sum(log(lik),na.rm=T)
  }
  return(llik)}




llike.GA = function(y,z,mu,alpha,p,rho){
  id0 = which(y==0)
  id1 = which(y!=0)

  u = pgam(y[id1],mu[id1],alpha,p[id1])
  u = ifelse(abs(u-.5) < .4999,u,sign(u-.5)*.4999+.5)
  llik = sum(-.5*(rho[id1]^2*qnorm(u)^2 - 2*qnorm(u)*z[id1]*rho[id1])/
               (1-rho[id1]^2),na.rm=T)
  llik = llik + sum(log(1-p[id1]) +
                      dgamma(y[id1],mu[id1]*alpha,alpha,log=T),na.rm=T)

  if(length(id0)>0){
    llik = llik + sum(pnorm(qnorm(p[id0]),mean=rho[id0]*z[id0],
                            sd=sqrt(1-rho[id0]^2),log.p=T),na.rm=T)
  }
  return(llik)}




llike.Beta = function(y,z,mu,alpha,p,rho){
  id0 = which(y==0)
  id1 = which(y!=0)

  u = pbet(y[id1],mu[id1],alpha,p[id1])
  u = ifelse(abs(u-.5) < .4999,u,sign(u-.5)*.4999+.5)
  llik = sum(-.5*(rho[id1]^2*qnorm(u)^2 - 2*qnorm(u)*z[id1]*rho[id1])/
               (1-rho[id1]^2),na.rm=T)
  llik = llik + sum(log(1-p[id1]) +
                      dbeta(y[id1],mu[id1]*alpha,(1-mu[id1])*alpha,log=T),na.rm=T)

  if(length(id0)>0){
    llik = llik + sum(pnorm(qnorm(p[id0]),mean=rho[id0]*z[id0],
                            sd=sqrt(1-rho[id0]^2),log.p=T),na.rm=T)
  }
  return(llik)}




beta.update = function(marginal,beta,y,x,z,mu.old,alpha,p,rho,B=200,h=200){

  linkfun <- meanlink_finder(marginal)
  llike <- llike_finder(marginal)

  n = nrow(x); d = ncol(x)
  it = if(is.vector(beta)){1} else {nrow(beta)}
  beta.old = if(is.vector(beta)){beta} else{beta[it,]}

  if(it <= B){
    cov.beta = diag(d)*.01
  } else {cov.beta = cov(beta[(it-h+1):it,]) + .0001*diag(d)}
  beta.star = beta.old + c(rnorm(d)%*%chol(cov.beta))

  mu.star = c(linkfun(x%*%beta.star))

  log.r = llike(y,z,mu.star,alpha,p,rho) +
    sum(dnorm(beta.star,0,100,log=T)) -
    llike(y,z,mu.old,alpha,p,rho) -
    sum(dnorm(beta.old,0,100,log=T))

  if((!is.na(log.r)) & (log.r > log(runif(1)))){
    return(list(beta=beta.star,mu=mu.star))
  } else{return(list(beta=beta.old,mu=mu.old))}
}



alphaNB.update = function(alpha,y,z,mu,p,rho,B=200,h=200){

  it = length(alpha); alpha.old = alpha[it]
  sd.alpha = ifelse(it <= B, .1, sd(log(alpha[(it-h+1):it])))
  alpha.star = exp(rnorm(1,log(alpha.old),sd.alpha))

  log.r = llike.NB(y,z,mu,alpha.star,p,rho) +
    dgamma(alpha.star,.001,.001,log=T) - log(alpha.old) -
    llike.NB(y,z,mu,alpha.old,p,rho) -
    dgamma(alpha.old,.001,.001,log=T) + log(alpha.star)

  if((!is.na(log.r)) & (log.r > log(runif(1)))){
    return(alpha.star)
  } else {return(alpha.old)}
}



alphaGA.update = function(alpha,y,z,mu,p,rho,B=200,h=200){

  it = length(alpha); alpha.old = alpha[it]
  sd.alpha = ifelse(it <= B, .1, sd(log(alpha[(it-h+1):it])))
  alpha.star = exp(rnorm(1,log(alpha.old),sd.alpha))

  log.r = llike.GA(y,z,mu,alpha.star,p,rho) +
    dgamma(alpha.star,.001,.001,log=T) - log(alpha.old) -
    llike.GA(y,z,mu,alpha.old,p,rho) -
    dgamma(alpha.old,.001,.001,log=T) + log(alpha.star)

  if((!is.na(log.r)) & (log.r > log(runif(1)))){
    return(alpha.star)
  } else {return(alpha.old)}
}



alphaBeta.update = function(alpha,y,z,mu,p,rho,B=200,h=200){

  it = length(alpha); alpha.old = alpha[it]
  sd.alpha = ifelse(it <= B, .1, sd(log(alpha[(it-h+1):it])))
  alpha.star = exp(rnorm(1,log(alpha.old),sd.alpha))

  log.r = llike.Beta(y,z,mu,alpha.star,p,rho) +
    dgamma(alpha.star,.01,.01,log=T) - log(alpha.old) -
    llike.Beta(y,z,mu,alpha.old,p,rho) -
    dgamma(alpha.old,.01,.01,log=T) + log(alpha.star)

  if((!is.na(log.r)) & (log.r > log(runif(1)))){
    return(alpha.star)
  } else {return(alpha.old)}
}



eta.update = function(marginal,eta,y,x,z,mu,alpha,p.old,rho,B=200,h=200){

  llike <- llike_finder(marginal)

  n = nrow(x); d = ncol(x)
  it = if(is.vector(eta)){1} else{nrow(eta)}
  eta.old = if(is.vector(eta)){eta} else{eta[it,]}

  if (!ZImarg(marginal)){
    return(list(eta=eta.old,p=p.old))
  } else {

    if(it <= B){
      cov.eta = diag(d)*.01
    } else{cov.eta = cov(eta[(it-h+1):it,]) + .001*diag(d)}
    eta.star = eta.old + c(rnorm(d)%*%chol(cov.eta))

    p.star = c(exp(x%*%eta.star)/(1 + exp(x%*%eta.star)))

    log.r = llike(y,z,mu,alpha,p.star,rho) +
      sum(dnorm(eta.star,0,100,log=T)) -
      llike(y,z,mu,alpha,p.old,rho) -
      sum(dnorm(eta.old,0,100,log=T))

    if((!is.na(log.r)) & (log.r > log(runif(1)))){
      return(list(eta=eta.star,p=p.star))
    } else{return(list(eta=eta.old,p=p.old))}
  }
}



zNB.update = function(y,z,mu,alpha,p,rho){

  uu = pnbinom(y,mu=mu,size=alpha)*(1-p) + p
  lu = pnbinom(y-1,mu=mu,size=alpha)*(1-p) + p*((y-1) >= 0)
  lu = ifelse(lu < uu, lu, max(uu - .0001,0))
  ub = qnorm(uu)
  lb = qnorm(lu)
  # lb = ifelse(lb<ub,lb,ub-.0001)
  z.draw = rtnorm(length(y),z*rho,sqrt(1-rho^2),lower=lb,upper=ub)
  return(z.draw)}



zGA.update = function(y,z,mu,alpha,p,rho){
  id0 = which(y==0)
  id1 = which(y!=0)
  z.draw = rep(0,length(y))

  u = pgam(y[id1],mu[id1],alpha,p[id1])
  u = ifelse(abs(u-.5) < .4999,u,sign(u-.5)*.4999+.5)
  z.draw[id1] = qnorm(u)

  if(length(id0)>0){
    z.draw[id0] = rtnorm(length(id0),mean=rho[id0]*z[id0],
                         sd=sqrt(1-rho[id0]^2),upper=c(qnorm(p[id0])))
  }
  return(z.draw)}



zBeta.update = function(y,z,mu,phi,p,rho){
  id0 = which(y==0)
  id1 = which(y!=0)
  z.draw = rep(0,length(y))

  u = pbet(y[id1],mu[id1],phi,p[id1])
  u = ifelse(abs(u-.5) < .4999,u,sign(u-.5)*.4999+.5)
  z.draw[id1] = qnorm(u)

  if(length(id0)>0){
    z.draw[id0] = rtnorm(length(id0),mean=rho[id0]*z[id0],
                         sd=sqrt(1-rho[id0]^2),upper=c(qnorm(p[id0])))
  }
  return(z.draw)}



tau.update = function(tau,z,x,rho.old,B=200,h=200){
  n = nrow(x); d = ncol(x)
  it = if(is.vector(tau)){1} else {nrow(tau)}
  tau.old = if(is.vector(tau)){tau} else{tau[it,]}

  if(it<=B){
    cov.tau = diag(d)*.01
  } else{cov.tau = cov(tau[(it-h+1):it,]) + .001*diag(d)}
  tau.star = tau.old + c(rnorm(d)%*%chol(cov.tau))

  rho.star = c((exp(x%*%tau.star)-1)/(exp(x%*%tau.star)+1))
  # rho.old = c((exp(x%*%tau.old)-1)/(exp(x%*%tau.old)+1))

  log.r = mvdens.z(z,rho.star) + sum(dnorm(tau.star,0,100,log=T)) -
    mvdens.z(z,rho.old) - sum(dnorm(tau.old,0,100,log=T))

  if((!is.na(log.r)) & (log.r > log(runif(1)))){
    return(list(tau=tau.star,rho=rho.star))
  } else{return(list(tau=tau.old,rho=rho.old))}
}

