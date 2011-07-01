\name{acdi}
\alias{acdi}
\title{Acceleration Deceleration}
\description{Acceleration Deceleration Technical Indicator}

\usage{acdi(Close, High = NULL, Low = NULL, Vol = NULL, plot = TRUE, ...)}

\arguments{\item{Close}{ VECTOR. Close price.}\item{High}{ VECTOR. High price.}\item{Low}{ VECTOR.
Low price.}\item{Vol}{ VECTOR. Asset traded Volume.}\item{plot}{ LOGICAL. If TRUE plot is
returned.}\item{...}{ Further arguments to or from other methods.}}

\details{Acceleration/-Deceleration Technical Indicator (AC) measures acceleration and deceleration
of the current driving force.  \cr 
This indicator will change direction before any changes in the driving force, which, it its turn,
will change its direction before the price.  \cr 
The nought line is basically the spot where the driving force is at balance with the acceleration. 
\cr 
If Acceleration-Deceleration is higher than nought, then it is usually easier for the acceleration
to continue the upward movement (and vice versa in cases when it is below nought). }


\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
