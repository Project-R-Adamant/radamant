\name{Appraisal}
\alias{Appraisal}
\alias{Appraisal.default}
\alias{Appraisal.Capm}
\title{Appraisal ratio}
\description{
\code{Appraisal}: Calculate Jensen index for a portfolio \cr
\code{Appraisal.Capm}: Get Jensen index from an object of class "Capm".
}
\usage{
Appraisal(PTF, ...)
\method{Appraisal}{default}(PTF, PTF_M, rf = NULL, rfr = 0, ...) 
\method{Appraisal}{Capm}(PTF, rfr = 0, ...)
}
\arguments{
\item{PTF}{ Input portfolio or an object of class "Capm"}
\item{PTF_M}{ Market/benchmark portfolio}
\item{rfr}{ risk free rate}
\item{rf}{ risk free asset}
\item{...}{ Further arguments to or from other methods}}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{Sharpe}},  \code{\link{Treynor}}, \code{\link{Jensen}}}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
