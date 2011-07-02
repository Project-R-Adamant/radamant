\name{treynor}
\alias{Treynor}
\alias{Treynor.default}
\alias{Treynor.Capm}
\title{Treynor index}
\description{
\code{Treynor}: Calculate Treynor index for a portfolio \cr
\code{Treynor.Capm}: Get Treynor index from an object of class "Capm"
}
\usage{
Treynor(PTF, ...)
\method{Treynor}{default}(PTF, PTF_M, rfr = 0, rf = NULL, ...) 
\method{Treynor}{Capm}(PTF, rfr = 0, ...)
}
\arguments{ 
\item{PTF}{ Input portfolio or an object of class "Capm"}
\item{PTF_M}{ Market/benchmark portfolio}
\item{rfr}{ risk free rate}
\item{rf}{ risk free asset}
\item{...}{ Further arguments to or from other methods}
}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{Jensen}}, \code{\link{Sharpe}}, \code{\link{Appraisal}} }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
