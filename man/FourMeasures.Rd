\name{FourMeasures}
\alias{FourMeasures}
\alias{FourMeasures.default}
\alias{FourMeasures.Capm}

\title{Four Measures indexes}

\description{Calculate the Four Measures indexes}

\usage{FourMeasures(PTF, PTF_M, rf = NULL, rfr = 0, ...) \cr
FourMeasures.Capm(X, rfr = 0, ...) }

\arguments{ 
\item{PTF}{ Input portfolio}
\item{PTF_M}{ Market/benchmark portfolio}
\item{X}{An object of class "Capm"}
\item{rfr}{ risk free rate}
\item{rf}{ risk free asset}
\item{...}{ }
}

\details{}

\value{Return a matrix containing the values for the following indexes: Sharpe, Treynor, Jensen and Appraisal}

\references{}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\note{}

\seealso{\code{\link{Sharpe}},  \code{\link{Treynor}}, \code{\link{Jensen}}, \code{\link{Appraisal}} }

\examples{}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
