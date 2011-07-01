\name{dpo}
\alias{dpo}
\title{Detrended price oscillator}
\description{Compute Detrended price oscillator (Technical Analysis)}
\usage{dpo(Close, lag = 5, plot = TRUE, ...) }
\arguments{\item{Close}{ VECTOR. Close price.}\item{lag}{ INTEGER. Number of lag
periods.}\item{plot}{ LOGICAL. If TRUE plot is returned.}\item{...}{ Further arguments to or from
other methods.}}

\details{The Detrended Price Oscillator removes the trend in prices by subtracting a moving average
of the price from the price. \cr
The Detrended Price shows cycles and overbought/oversold conditions. Note that the calculation
shifts the results (shift = term / 2 + 1) periods, so the last shift periods will be zero.}


\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
