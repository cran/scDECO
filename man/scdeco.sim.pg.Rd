% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scdeco.sim.pg.R
\name{scdeco.sim.pg}
\alias{scdeco.sim.pg}
\title{Simulating from ZENCO Model}
\usage{
scdeco.sim.pg(
  N,
  b0,
  b1,
  phi1,
  phi2,
  mu1,
  mu2,
  tau0,
  tau1,
  mu3,
  phi3,
  tau2 = NULL,
  tau3 = NULL,
  xc = NULL
)
}
\arguments{
\item{N}{size of sample to be generated}

\item{b0}{intercept of zinf parameter}

\item{b1}{slope of zinf parameter}

\item{phi1}{over-dispersion parameter of first marginal}

\item{phi2}{over-dispersion parameter of second marginal}

\item{mu1}{mean parameter of first marginal}

\item{mu2}{mean parameter of second marginal}

\item{tau0}{intercept of correlation}

\item{tau1}{slope of of correlation}

\item{mu3}{mean parameter of covariate vector}

\item{phi3}{over-dispersion parameter of covariate vector}

\item{tau2}{(optional) correlation coefficient on optional xc covariate vector}

\item{tau3}{(optional) correlation coefficient on interaction between x3 and xc}

\item{xc}{(optional) secondary covariate to be regressed}
}
\value{
a matrix with expressions as first two columns and covariates as remaining columns
}
\description{
Simulating from ZENCO Model
}
\examples{
phi1_use <- 4
phi2_use <- 4
phi3_use <- 1/6
mu1_use <- 15
mu2_use <- 15
mu3_use <- 7
b0_use <- 0.6882
b1_use <- -0.2995
tau0_use <- 0.07
tau1_use <- 0.05

simdat <- scdeco.sim.pg(N=1000, b0=b0_use, b1=b1_use,
                        phi1=phi1_use, phi2=phi2_use, phi3=phi3_use,
                        mu1=mu1_use, mu2=mu2_use, mu3=mu3_use,
                        tau0=tau0_use, tau1=tau1_use)
simdat[1:10,]
}
