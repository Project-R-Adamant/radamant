\name{Jensen}
\alias{Jensen}
\alias{Jensen.default}
\alias{Jensen.Capm}

\title{Jensen index}

\description{
\code{Jensen}: Calculate Jensen index for a portfolio \cr
\code{Jensen.Capm}: Get Jensen index from an object of class "Capm".
}

\usage{
Jensen(x, ...)

\method{Jensen}{default}(PTF, PTF_M, rf = NULL, rfr = 0, ...) 

\method{Jensen}{Capm}(X, rfr = 0, ...)

}

\arguments{ 
\item{PTF}{ Input portfolio}
\item{PTF_M}{ Market/benchmark portfolio}
\item{X}{An object of class "Capm"}
\item{rfr}{ risk free rate}
\item{rf}{ risk free asset}
\item{...}{ }}

\details{}

\value{}

\references{}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\note{}

\seealso{\code{\link{Sharpe}},  \code{\link{Treynor}}, \code{\link{Appraisal}}}

\examples{}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
