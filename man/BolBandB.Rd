\name{BolBandB}

\alias{BolBandB}

\title{Bollinger Bands Bandwidth}

\description{Compute Bollinger Bands Bandwidth (Technical analysis)}
\usage{BolBandB(Close, High, Low, fact=2, win.size=5, plot=FALSE, ...)}
\arguments{\item{Close}{ VECTOR. Close price.}\item{High}{ VECTOR. High price.}\item{Low}{ VECTOR.
Low price.}\item{fact}{ fact}\item{win.size}{ win.size}\item{plot}{ LOGICAL. If TRUE plot is
returned.}\item{...}{ Further arguments to or from other methods.}}
\details{Bollinger Bandwidth is an important indicator derived from John Bollingers original
Bollinger Bands indicator.  \cr Bandwidth is a relative measure of the width of the Bollinger Bands.
\cr The Bollinger Bandwidth equation can be constructed as: \cr Volatility is high when the
Bollinger Bands are farther apart and low when the Bollinger Bands are closer together. Based on the
assumption that price (and volatility) generally operates in cycles - periods of low volatility
inevitably followed by periods of high volatility and so on and so forth - traders can learn to take
advantage of these cycles. \cr The Bollinger Band Width indicator is the distance between the upper
and lower Bollinger Bands. It is a measure of volatility. The Band Width value is higher when
volatility is high, and lower when volatility is low. High Band Width values indicate that the
current trend may be about to end. Low Band Width values indicate that a new trend may be about to
start. \cr Volatility is high when the Bollinger Bandwidth moves to its high levels and volatility.
It is low when Bandwidth moves closer to zero. As a rule, a stocks price and volatility moves in
cycles - the price moves up and down and periods of low volatility are replaced by periods of high
volatility. In general, periods of high volatility can be noted during down-trends and corrections
downward. Periods of low volatility can be observed during up-trends and recoveries. With this
knowledge a trader can build a trading system that generates signals based on the volatility
technical analysis. \cr The Bollinger Bandwidth can be used to identify the Squeeze - when the
Bandwidth is at its lowest low value within n-periods. Bollinger states that Squeeze could occur
before a trend reversal, like the calm before the storm. In this case, trading Buy/Sell signals can
be generated on the price breakouts following the Squeeze.}
\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
