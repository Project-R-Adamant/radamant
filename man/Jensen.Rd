\name{jensen}
\alias{Jensen}
\alias{Jensen.default}
\alias{Jensen.Capm}
\title{Jensen index}
\description{
\code{Jensen}: Calculate Jensen index for a portfolio.
\code{Jensen.Capm}: Get Jensen index from an object of class "Capm".
}
\usage{
Jensen(PTF, ...)
\method{Jensen}{default}(PTF, PTF_M, rf = NULL, rfr = 0, ...) 
\method{Jensen}{Capm}(PTF, rfr = 0, ...)
}
\arguments{ 
\item{PTF}{ Input portfolio or an object of class "Capm"}
\item{PTF_M}{ Market/benchmark portfolio}
\item{rfr}{ risk free rate}
\item{rf}{ risk free asset}
\item{...}{ Further arguments to or from other methods}
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{Sharpe}},  \code{\link{Treynor}}, \code{\link{Appraisal}}}


