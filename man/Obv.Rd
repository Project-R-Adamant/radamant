\name{obv}
\alias{Obv}
\title{On Balance Volume oscillator}
\description{Compute On Balance Volume oscillator (Technical Analysis)}
\usage{Obv(Close, Volume) }
\arguments{\item{Close}{ VECTOR. Close price.}\item{Volume}{ VECTOR. Asset traded Volume.}}

\details{
The On Balance Volume (OBV) is a cumulative total of the up and down volume. When the
close is higher than the previous close, the volume is added to the running total, and when the
close is lower than the previous close, the volume is subtracted from the running total. \cr
To interpret the OBV, look for the OBV to move with the price or precede price moves. If
the price moves before the OBV, then it is a non-confirmed move. A series of rising peaks, or
falling troughs, in the OBV indicates a strong trend. If the OBV is flat, then the market is not
trending.
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
