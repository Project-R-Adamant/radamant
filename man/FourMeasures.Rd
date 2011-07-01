\name{FourMeasures}

\alias{FourMeasures}
\alias{FourMeasures.default}
\alias{FourMeasures.Capm}

\title{Four Measures indexes}

\description{Calculate the Four Measures indexes}

\usage{
FourMeasures(PTF, ...)
\method{FourMeasures}{default}(PTF, PTF_M, rf = NULL, rfr = 0, ...) 
\method{FourMeasures}{Capm}(PTF, rfr = 0, ...) 
}

\arguments{ 
\item{PTF}{ Input portfolio or an object of class "Capm"}
\item{PTF_M}{ Market/benchmark portfolio}
\item{rfr}{ risk free rate}
\item{rf}{ risk free asset}
\item{...}{ Further arguments to or from other methods}
}

\value{Return a matrix containing the values for the following indexes: Sharpe, Treynor, Jensen and
Appraisal}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{Sharpe}},  \code{\link{Treynor}}, \code{\link{Jensen}}, \code{\link{Appraisal}}
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
