\name{Abi}

\alias{Abi}

\title{Absolute Breath Index - ABI}

\description{Compute Absolute Breath Index (Technical Analysis)}

\usage{Abi(X, lag = 5, plot=FALSE, ...)}

\arguments{ 
\item{X}{ Input numerical series}
\item{lag}{ Number of lags}
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ Further arguments to or from other methods}
}

\details{The absolute breadth index (ABI) measures the amount of volatile on the New York Stock
Exchange. The indicator is truly a momentum indicator since it only tracks the movement of the
issues on the exchange and not their direction. The overall purpose of the indicator is to determine
when the market is most volatile, thus providing the greatest opportunity for larger profits due to
increased volatility. \cr Calculation: \cr \deqn{ABI = |Advancing \ Issues - Declining \ issues|}}

\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
