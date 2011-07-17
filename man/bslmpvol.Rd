\name{bslmpvol}
\alias{BS.ImpVol}
\title{Black & Scholes Implied volatility}
\description{Calculate Black & Scholes Implied volatility}

\usage{
BS.ImpVol(P, under, strike, rfr, sigma, maty, 
yield, 
calc.type =c("standard", "lognorm", "gammarec"), 
opt.type = c("call", "put"),
interval = c(-20, 20)) 
}

\arguments{
\item{P}{ Observed Price; single numeric }
\item{under}{ Underlying asset price.}
\item{strike}{ Strike/Exercise price.}
\item{rfr}{ Risk free rate (continuos)}
\item{sigma}{ Assets standard deviation - annualised volatility.}
\item{maty}{ Period of maturity.}
\item{yield}{ Dividend yield (continuos)}
\item{calc.type}{ Calculation type.}
\item{opt.type}{ Type of option (Default="call").}
\item{interval}{ calculation interval applied to the function uniroot (\code{\link{uniroot}})}
}

\value{
Matrix of Px1 dimensions with Implied volatility values. One row for each value of P.
}
\author{RAdamant Development Team \email{team@r-adamant.org}}

\seealso{ \code{\link{BS.greeks}}, \code{\link{uniroot}}}

\examples{
	# Set BS paramaters
	under<- 100
	strike <- 95
	rfr<- 0.08
	sigma <- 0.2
	maty<-0.5
	yield<- 0.03
	calc.type<-"lognorm"
	opt.type<-"call"
	
	# calculate implied volatility for single oberved price
	P = 11
	imp = BS.ImpVol(P, under, strike, rfr, sigma, maty, yield)
	imp
	# calculate implied volatility for multiple oberved prices
	P = seq(9, 11, by=0.1)
	imp = BS.ImpVol(P, under, strike, rfr, sigma, maty, yield)
	imp
}
