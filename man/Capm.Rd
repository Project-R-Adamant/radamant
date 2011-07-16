\name{capm}

\alias{Capm.default}
\alias{Capm}

\title{Capm - default method}

\description{Default method for CAPM}

\usage{
Capm(PTF, ...)
\method{Capm}{default}(PTF, PTF_M, rf = NULL, rfr = NULL, ...) 
}

\arguments{ 
\item{PTF}{ Matrix of returns, one series for each asset in the portfolio.}
\item{PTF_M}{ Vector of returns for the market portfolio}
\item{rf}{ Vector. Risk free asset returns}
\item{rfr}{ Numeric. Risk free rate}
\item{...}{ Further arguments to or from other methods}}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{

	# load example dataset
	data(ex_ptf)
	# Generate a random return risk free asset
	rf = rnorm(NROW(ex_ptf), mean = 0.05, sd = 0.01)
	# Calculate CAPM
	Capm(PTF = ex_ptf[,-1], PTF_M = ex_ptf[,1], rf)

\dontrun{
	## Example with real time series 
	ACME = get.fs("APKT",SName = "Acme Packet", from=as.Date("2010-01-01"))
	ABTL = get.fs("ABTL", SName = "Autobytel", from=as.Date("2010-01-01"))
	CNAF = get.fs("CNAF", from=as.Date("2010-01-01"))
	BIIB = get.fs("BIIB", SName = "Biogen", from=as.Date("2010-01-01"))
	SONY = get.fs("SNE", SName = "Sony", from=as.Date("2010-01-01"))
	ENI = get.fs("E", SName = "Eni", from=as.Date("2010-01-01"))
	ptf = combine.fs(ACME, ABTL, CNAF, BIIB, SONY, ENI);
	head(ptf)

	# Load a Benchmark Portfolio Index
	NASDAQ = get.fs("^IXIC", SName = "NASDAQ", from=as.Date("2010-01-01"));

	R_ptf = Ret(ptf, na.rm = TRUE);
	# Return of the Benchmark portfolio (NASDAQ index)
	R_NASDAQ = Ret(NASDAQ, na.rm = TRUE)

	# Generate a random return risk free asset
	rf = rnorm(NROW(R_ptf), mean = 0.05, sd = 0.01)
	Capm(R_ptf, R_NASDAQ, rf)
}
}



