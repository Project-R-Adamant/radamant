\name{VecAr}

\alias{VecAr}
\alias{VecAr.default}

\title{Vector Autoregressive model}

\description{Estimate Vector Autoregressive model}

\usage{
VecAr(x, ...)
\method{VecAr}{default}(Y, ar.lags = 1:2, type = c("const", "trend", "constrend", "none"),
exog = NULL) 
}

\arguments{ 
\item{Y}{ Input matrix of time series. \cr
	N.B. The first column is taken as dependent variable }
\item{ar.lags}{ Number (or vector) of lags for the AR components }
\item{type}{ Type of deterministic regressor(s) to be included in the model}
\item{exog}{ matrix of exogenous variables (Default = NULL)}
}

\details{}

\value{An object list of class "VecAr". The list contains the following elements: \cr
	- Results of the estimation ("lm" object) \cr
	- Nunmber of Observations \cr
	- Number of Variables \cr
	- Number of Parameters \cr
	- LogLikelihood value \cr
	- AIC information criteria \cr
	- BIC information criteria
	}

\references{}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\note{}

\seealso{ \code{\link{Strvar}}, \code{\link{fitted.VecAr}}  }

\examples{}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
