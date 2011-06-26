\name{Treynor}
\alias{Treynor}
\alias{Treynor.default}
\alias{Treynor.Capm}

\title{Treynor index}

\description{
\code{Treynor}: Calculate Treynor index for a portfolio \cr
\code{Treynor.Capm}: Get Treynor index from an object of class "Capm"
}

\usage{
Treynor(x, ...)

\method{Treynor}{default}(PTF, PTF_M, rfr = 0, rf = NULL, ...) 

\method{Treynor}{Capm}(X, rfr = 0, ...)

}

\arguments{ 
\item{PTF}{ Input portfolio}
\item{PTF_M}{ Market/benchmark portfolio}
\item{rfr}{ risk free rate}
\item{rf}{ risk free asset}
\item{X}{An object of class "Capm"}
\item{...}{ }
}

\details{}

\value{}

\references{}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\note{}

\seealso{\code{\link{Jensen}}, \code{\link{Sharpe}}, \code{\link{Appraisal}} }

\examples{}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
