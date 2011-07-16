\name{garch}

\alias{Garch}
\alias{Garch.default}

\title{Garch}

\description{Estimate Garch models}

\usage{
Garch(x, ...)
\method{Garch}{default}(x, Y=NULL, order=c(alpha=1,beta=1), phi=0, delta=0,
type=c("garch","mgarch","tgarch","egarch"), prob=c("norm","ged","t"), ...)
}

\arguments{ 
\item{x}{ Univariate time series, usually returns}
\item{Y}{ Exogenous regressors for the Mean Equation}
\item{order}{ Garch order}
\item{type}{ Garch type.}
\item{prob}{ Probability density for the innovations.}
\item{phi}{ Phi pars}
\item{delta}{ Delta pars}
\item{...}{ Further arguments to or from other methods}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}


