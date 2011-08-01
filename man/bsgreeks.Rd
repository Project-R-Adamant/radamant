\name{bsgreeks}
\alias{BS.greeks}

\title{Black & Scholes greeks}

\description{Calculate analytically Black & Scholes greeks}

\usage{BS.greeks(X = NULL, ...) }

\arguments{ 
\item{X}{ An object of class "BS.price"}
\item{...}{ Further arguments to or from other methods - parameters accepted by the function 
 \code{\link{BS.price}}}
}

\value{
A matrix containing the values for calculated greeks:
\item{Delta}{}
\item{Vega}{}
\item{Theta}{}
\item{Rho}{}
\item{Lambda}{}
\item{Gamma}{}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\seealso{ \code{\link{BS.price}}, \code{\link{BS.moments}} }

\examples{
# Set BS paramaters
under = 105
strike = 95
rfr = 0.08
sigma = 0.2
maty = 0.5
yield = 0.03

# calculate BS price for a call option 
# assuming normal distribution of prices
bs1 = BS.price(under
		, strike
		, rfr
		, sigma
		, maty
		, yield
		, calc.type = "standard"
		, opt.type = "call") 

# assuming gamma-reciprocal distribution of prices
bs2 = BS.price(under
		, strike
		, rfr
		, sigma
		, maty
		, yield
		, calc.type = "gammarec"
		, opt.type = "call") 

# calculate greeks for object bs1 of class "BS.price"
BS.greeks(bs1)
class(bs1)

# ... or alternatively passing the same BS paramaters
#  used for price calculation the results are the same
BS.greeks(under=under
	, strike=strike
	, rfr=rfr
	, sigma=sigma
	, maty=maty
	, yield=yield
	, opt.type = "call"
	, calc.type = "standard" )

# Same examples as above for different calculation type
BS.greeks(bs2)
class(bs2)
BS.greeks(under=under
	, strike=strike
	, rfr=rfr
	, sigma=sigma
	, maty=maty
	, yield=yield
	, opt.type = "call"
	, calc.type = "gammarec")
}


