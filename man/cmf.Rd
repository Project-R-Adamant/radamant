\name{cmf}
\alias{cmf}
\title{Chaikin Money Flow}
\description{Compute Chaikin Money Flow (Technical Analysis)}
\usage{cmf(Close, Low, High, Volume, plot = FALSE, ...) }

\arguments{ \item{Close}{ VECTOR. Close price.}\item{Low}{ VECTOR. Low price.}\item{High}{ VECTOR.
High
price.}\item{Volume}{ Volume}\item{plot}{ LOGICAL. If TRUE plot is returned.}\item{...}{ Further
arguments to or from other
methods}}

\details{Chaikin Money Flow measures the amount of Money Flow Volume over a specific period. \cr 
Money Flow Volume forms the basis for the Accumulation Distribution Line. Chaikin Money Flow sums
Money Flow Volume for a specific look-back period, typically 20 or 21 days. The resulting indicator
fluctuates above/below the zero line just like an oscillator. Chartists weigh the balance of buying
or selling pressure with the absolute level of Chaikin Money Flow. Chartists can also look for
crosses above or below the zero line to identify changes on money flow. The Chaikin Money Flow
compares the total volume over the last n time periods to the total of volume times the Closing
Location Value (CLV) over the last n time periods. The CLV calculates where the issue closes within
its trading range. When the Chaikin Money Flow is above 0.25 it is a bullish signal, when it is
below -0.25, it is a bearish signal. If the Chaikin Money Flow remains below zero while the price is
rising, it indicates a probable reversal.}


\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
