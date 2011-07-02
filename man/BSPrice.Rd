\name{bsprice}

\alias{BS.price.default}
\alias{BS.price}

\title{Black & Scholes price generic}

\description{Generic method for Black & Scholes price}

\usage{
BS.price(under, ...)
\method{BS.price}{default}(under, strike, rfr, sigma, maty, yield, calc.type =c("standard",
"lognorm","gammarec"), opt.type = c("call", "put"), ...) 
}

\arguments{ 
\item{under}{ Underlying asset price.}
\item{strike}{ Strike price.}
\item{rfr}{ Risk free rate.}
\item{sigma}{ Volatility.}
\item{maty}{ Maturity.}
\item{yield}{ Yield}
\item{calc.type}{ Calculation type.}
\item{opt.type}{ Option type.}
\item{...}{ Further arguments to or from other methods}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
