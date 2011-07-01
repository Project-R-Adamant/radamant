\name{VecAr}

\alias{VecAr}
\alias{VecAr.default}

\title{Vector Autoregressive model}

\description{Estimate Vector Autoregressive model}

\usage{
VecAr(X, ...)

\method{VecAr}{default}(X, ar.lags = 1:2, 
type = c("const", "trend", "constrend", "none"),
exog = NULL, ...) 
}

\arguments{ 
\item{X}{ Input matrix of time series. \cr
	N.B. The first column is taken as dependent variable }
\item{ar.lags}{ Number (or vector) of lags for the AR components }
\item{type}{ Type of deterministic regressor(s) to be included in the model}
\item{exog}{ matrix of exogenous variables (Default = NULL)}
\item{...}{ Further arguments to or from other methods}
}

\value{An object list of class "VecAr". The list contains the following elements: \cr
	- Results of the estimation ("lm" object) \cr
	- Nunmber of Observations \cr
	- Number of Variables \cr
	- Number of Parameters \cr
	- LogLikelihood value \cr
	- AIC information criteria \cr
	- BIC information criteria
	}

\author{ RAdamant Development Team \email{team@r-adamant.org} }
	
\seealso{ \code{\link{Strvar.VecAr}}, \code{\link{fitted.VecAr}}  }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
