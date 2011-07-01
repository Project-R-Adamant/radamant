\name{cmof}
\alias{cmof}
\title{Chande Momentum Oscillator}
\description{Compute Chande Momentum Oscillator (Technical Analysis)}
\usage{cmof(X, lag = 5, plot = FALSE, ...) }
\arguments{\item{X}{ X}\item{lag}{ INTEGER. Number of lag periods.}\item{plot}{ LOGICAL. If TRUE
plot is returned.}\item{...}{ Further arguments to or from other methods.}}

\details{The Chande Momentum Oscillator is a modified RSI. Where the RSI divides the upward movement
by the net movement (up / (up + down)), the CMO divides the total movement by the net movement ((up
- down) / (up + down)). There are several ways to interpret the CMO. Values over 50 indicate
overbought conditions, while values under -50 indicate oversold conditions. High CMO values indicate
strong trends. When the CMO crosses above a moving average of the CMO, it is a buy signal, crossing
down is a sell signal.}


\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
