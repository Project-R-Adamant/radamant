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
Appraisal(x, ...)

\method{Appraisal}{default}(PTF, PTF_M, rf = NULL, rfr = 0, ...) 

\method{Appraisal}{Capm}(X, rfr = 0, ...)

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

\seealso{\code{\link{Sharpe}},  \code{\link{Treynor}}, \code{\link{Jensen}}}

\examples{}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
