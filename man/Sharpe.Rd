\name{sharpe}
\alias{Sharpe}
\alias{Sharpe.default}
\alias{Sharpe.Capm}
\title{Sharpe index}

\description{
\code{Sharpe}: Calculate Sharpe index for a portfolio. \cr
\code{Sharpe.Capm}: Get Sharpe index from an object of class. "Capm"
}

\usage{
Sharpe(PTF, ...)
\method{Sharpe}{default}(PTF, rfr = 0, ...) 
\method{Sharpe}{Capm}(PTF, rfr = 0, ...) 
}

\arguments{ 

\item{PTF}{ Input portfolio or an object of class "Capm"}

\item{rfr}{ risk free rate}

\item{...}{ Further arguments to or from other methods}
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{Treynor}}, \code{\link{Jensen}}, \code{\link{Appraisal}}}


