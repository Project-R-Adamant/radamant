\name{ptfopt}

\alias{PtfOpt}
\alias{print.PtfOpt}

\title{Mean-Variance optimum portfolio}

\description{Calculate mean-variance efficient portfolio}

\usage{
PtfOpt(ret = NULL
	, ptf = NULL
	, mi = NULL
	, SIGMA = NULL
	, volatility = TRUE
	, ...
	)
    
\method{print}{PtfOpt}(x, ...) 
}

\arguments{
\item{ret}{ Vector containing averge return for each asset}
\item{ptf}{ Matrix containing one or more series of prices, one time series for each asset}
\item{mi}{ Target return for the portfolio}
\item{SIGMA}{ Sample covariance matrix}
\item{volatility}{ Logical. If TRUE volatility is returned, else the variance is computed.}
\item{x}{ An object of class "PrfOpt".}
\item{...}{ Further arguments to or from other methods}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\seealso{\code{\link{PtfFront}}, \code{\link{PtfUtility}} }

\examples{

# Calculate weights from a series of prices
data(EuStockMarkets);
PtfOpt(ptf = EuStockMarkets);

# simulate efficient frontier
PtfFront(PTF = EuStockMarkets
	, n_sim=100
	, col="yellow"
	);
PtfFront(PTF = EuStockMarkets
	, n_sim=30
	, col="green"
	);

# calculate weights from a vector of returns R and matrix SIGMA
R = c(A = 0.021, B = 0.09);
SIGMA = matrix(c(0.101^2, 0.005,  0.005, 0.208^2)
               , nrow = 2
               , ncol = 2
               );
# set target returns to be 0.05
PtfOpt(ret = R
	, ptf = NULL
	, SIGMA = SIGMA
	, mi = c(0.05)
	);

# set two target returns: 0.05 and 0.07
PtfOpt(ret = R
	, ptf = NULL
	, SIGMA = SIGMA
	, mi = c(0.05, 0.07)
	);
# simulate efficient frontier
PtfFront(ret = R
	 , ptf = NULL
	 , SIGMA = SIGMA
	 , n_sim=100
	 , col="yellow"
	 );


## Example with real time series
fromDt = as.Date("2010-01-01");
ACME = get.fs("APKT",SName = "Acme Packet", from = fromDt);
ABTL = get.fs("ABTL", SName = "Autobytel", from = fromDt);
CNAF = get.fs("CNAF", from = fromDt);
BIIB = get.fs("BIIB", SName = "Biogen", from = fromDt);
SONY = get.fs("SNE", SName = "Sony", from = fromDt);
ENI = get.fs("E", SName = "Eni", from = fromDt);
ptf = combine.fs(ACME, ABTL, CNAF, BIIB, SONY, ENI);
head(ptf);

# Compute Minimum Variance portfolio
PtfOpt(ptf = ptf);

}



