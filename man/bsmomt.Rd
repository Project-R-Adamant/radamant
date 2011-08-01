\name{bsmomt}

\alias{BS.moments}

\title{Black & Scholes moments}

\description{Calculate first four moments for Black & Scholes}

\usage{BS.moments(BS = NULL, under, rfr, sigma, yield, maty) }

\arguments{
\item{BS}{ An object of class "BS.price"}
\item{under}{ Underlying asset price.}
\item{rfr}{ Risk free rate (continuos)}
\item{sigma}{ Assets standard deviation - annualised volatility.}
\item{yield}{ Dividend yield (continuos)}
\item{maty}{ Period of maturity.}
}

\value{
A matrix containing the four moments (one for each row):
\item{Mom_1}{}
\item{Mom_2}{}
\item{Mean}{}
\item{Var}{}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
    # Set BS paramaters
    under = 105
    strike = 95
    rfr = 0.08
    sigma = 0.2
    maty = 0.5
    yield = 0.03
    # calculate BS price
    bs = BS.price(under, strike, rfr, sigma, maty, yield)
    # calculate moments for object bs of class "BS.price"
    BS.moments(bs)
    # ... or alternatively passing the same BS paramaters used for price 
    # calculation, the results are the same
    BS.moments(NULL, under, rfr, sigma, yield, maty)

}


