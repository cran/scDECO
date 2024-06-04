#' ZENCO fitting function when secondary covariate is provided
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
#' @import rjags
#'
#' @keywords internal
#' @return MCMC samples that have been adapted, burned, and thinned
#'
scdeco.pg.xc <- function(dat, b0, b1, adapt_iter=10000, update_iter=5000, coda_iter=50000, coda_thin=20, coda_burnin=10000){
  ### JAGS ###
  IndZ.string <-"
  model{
  for (i in 1:N){
  p[i, 1] ~ dbern(exp(b0+b1*mu[1])/(1+exp(b0+b1*mu[1])))
  p[i, 2] ~ dbern(exp(b0+b1*mu[2])/(1+exp(b0+b1*mu[2])))
  p3[i] ~ dbern(exp(b0+b1*mu3)/(1+exp(b0+b1*mu3)))

  x[i] ~ dpois(lmd3[i])
  h[i] ~ dgamma(inverphi3, inverphi3)
  lmd3[i] <- ifelse(p3[i]==0, h[i]*mu3, 1e-10)


  rho_tmp[i] <- ifelse(p3[i]==0, ((exp(tau0 +tau1*x[i] +tau2*c[i] + tau3*x[i]*c[i])-1)/(exp(tau0+tau1*x[i]+tau2*c[i] + tau3*x[i]*c[i])+1)), (exp(tau0+tau1*mu3+tau2*c[i] + tau3*x[i]*c[i])-1)/(exp(tau0+tau1*mu3+tau2*c[i] + tau3*x[i]*c[i])+1))
  rho[i] <- ifelse(abs(rho_tmp[i])==1, rho_tmp[i]/(1+10^-9), rho_tmp[i])
  sigma[1,1,i] <- 1/(1-rho[i]*rho[i])
  sigma[2,2,i] <- 1/(1-rho[i]*rho[i])
  sigma[1,2,i] <- -rho[i]/(1-rho[i]*rho[i])
  sigma[2,1,i] <- -rho[i]/(1-rho[i]*rho[i])
  Z[i,1:2] ~ dmnorm(rep(0,2),sigma[,,i])

  for (j in 1:2){
  Y[i,j] ~ dpois(lambda[i,j])
  lambda[i,j] <- ifelse(p[i,j]==0, qgamma(pnorm(Z[i,j],0,1),inverphi[j],inverphi[j])*mu[j], 1e-10)
  }
  }
  tau0 ~ dnorm(0, 4/N)
  tau1 ~ dnorm(0, 4/N)
  tau2 ~ dnorm(0, 4/N)
  tau3 ~ dnorm(0, 4/N)
  mu3 ~ dlnorm(0, 1)
  inverphi3 ~ dgamma(1, 0.01)
  for (i in 1:2){
  mu[i] ~ dlnorm(0, 1)
  inverphi[i] ~ dgamma(1, 0.01)
  }
  }
  "
  y.indz <- cbind(dat[,1], dat[,2])
  x3.indz <- dat[,3]
  c.indz <- dat[,4]
  IndZ.spec <- textConnection(IndZ.string)
  jags_data = list(N=length(dat[,1]),Y=y.indz, x=x3.indz, c=c.indz, b0=b0, b1=b1)
  phi.indz <- c((var(y.indz[,1])/mean(y.indz[,1])-1)/mean(y.indz[,1]), (var(y.indz[,2])/mean(y.indz[,2])-1)/mean(y.indz[,2]))
  phi.indz[1]<-ifelse(phi.indz[1]<0, 1e-10, phi.indz[1])
  phi.indz[2]<-ifelse(phi.indz[2]<0, 1e-10, phi.indz[2])
  jags_inits = list(mu=c(mean(y.indz[,1]), mean(y.indz[,2])), inverphi=1/phi.indz, mu3=mean(x3.indz), inverphi3=1/((var(x3.indz)/mean(x3.indz)-1)/mean(x3.indz)),tau0=0, tau1=0, tau2=0, tau3=0)
  jags_model_IndZ = jags.model(IndZ.spec, data=jags_data, n.adapt=adapt_iter, inits=jags_inits, n.chains=3)
  update(jags_model_IndZ, update_iter)
  samps.coda.IndZ <- coda.samples(jags_model_IndZ, c('mu', 'inverphi', 'tau0', 'tau1', 'tau2', 'tau3', 'mu3', 'inverphi3'), n.iter = coda_iter, thin=coda_thin, n.burnin=coda_burnin)

  # summary(samps.coda.IndZ)
  samps.IndZ<- summary(samps.coda.IndZ)
  return(samps.IndZ)
}
