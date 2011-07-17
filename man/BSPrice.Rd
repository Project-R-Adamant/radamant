\name{bsprice}

\alias{BS.price.default}
\alias{BS.price}
\alias{print.BS.price}

\title{Black & Scholes price generic}

\description{Generic method for Black & Scholes price}

\usage{
BS.price(under, ...)
\method{BS.price}{default}(under
	, strike
	, rfr
	, sigma
	, maty
	, yield
	, calc.type =c("standard","lognorm","gammarec")
	, opt.type = c("call", "put")
	, ...) 
\method{print}{BS.price}(x, mod, ...)	
}

\arguments{ 
\item{under}{ Underlying asset price.}
\item{strike}{ Strike/Exercise price.}
\item{rfr}{ Risk free rate (continuos)}
\item{sigma}{ Assets standard deviation - annualised volatility.}
\item{maty}{ Period of maturity.}
\item{yield}{ Dividend yield (continuos)}
\item{calc.type}{ Calculation type.}
\item{opt.type}{ Type of option (Default="call").}
\item{x}{ An object of class "BS.price".}
\item{mod}{ Control object for print method.}
\item{...}{ Further arguments to or from other methods.}
}

\value{
An object of class "BS.price" containing:
\item{BS Price}{}
\item{Factor d1}{}
\item{Factor d2}{}
}

\details{
The parameter "calc.type" allows to change the Black & Scholes calculation according to different distributional
assumptions. 
\itemize{
\item{\code{standard}: Log asset price normally distributed}
\item{\code{lognorm}: Log asset price log-normally distributed}
\item{\code{gammarec}: Log asset price Gamma-Reciprocal distributed}
}
}     


\author{RAdamant Development Team \email{team@r-adamant.org}}

\seealso{ \code{\link{BS.greeks}}, \code{\link{BS.moments}} }

\examples{
# Set BS paramaters
under = 100
strike = 95
rfr = 0.08
sigma = 0.2
maty = 0.5
yield = 0.03

# calculate BS price for a call option assuming normal distribution of prices
bs1 = BS.price(under, strike, rfr, sigma, maty, yield, calc.type = "standard" , opt.type = "call") 
bs1
# same example assuming gamma-reciprocal distribution of prices
bs2 = BS.price(under, strike, rfr, sigma, maty, yield, calc.type = "gammarec" , opt.type = "call") 
bs2
}


