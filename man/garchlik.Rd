\name{garchlik}

\alias{like.garch}
\alias{like.tgarch}
\alias{like.egarch}
\alias{like.mgarch}

\title{GARCH likelihood functions}

\description{Calculate likelihood for Garch, TGarch, EGarch and MGarch models}

\usage{
like.mgarch(theta, x, Y, order, k, prob=c("norm","ged","t")) 
like.garch(theta, ee, x, Y, order, k, prob = c("norm","ged", "t")) 
like.tgarch(theta, ee, x, Y, order, k, prob = c("norm","ged", "t")) 
like.egarch(theta, ee, x, Y, order, k, prob = c("norm","ged", "t")) 
}

\arguments{
\item{theta}{ Vector of paramaters.}
\item{ee}{ Vector of innovations.}
\item{x}{ Original series of returns.}
\item{k}{ Number of mean equation regressors.}
\item{Y}{ Matrix of exogenous variables used for the mean equation.}
\item{order}{ Model parameter order}
\item{prob}{ Innovations probability density: "norm", "ged", "t". (Default: "norm")}
}

\value{Likelihood value}

\details{Those functions are called inside the main Garch function in order to obtain numerical optimisation of the input parameters. \cr
The input parameter of the functions are calculated directly inside the Garch function (see \code{\link{Garch}})
}

\author{RAdamant Development Team \email{team@r-adamant.org}}


