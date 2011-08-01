\name{styles}

\alias{Styles}
\alias{Multi.Styles}

\title{Styles analysis (portfolio)}

\description{Perform Style analysis for single and multiple time periods}

\usage{
Styles(FUND
	, IND
	, W
	, lower = NULL
	, upper = NULL
	, ...)

Multi.Styles(FUND
	, IND
	, W
	, n_clust = 5
	, lower = NULL
	, upper = NULL
	, ...) 
}

\arguments{
\item{FUND}{ Vector. Benchmark investment fund }
\item{IND}{ Matrix of indices (returns)}
\item{W}{ Initial weghts to be assigned to the indices}
\item{n_clust}{ Number of time periods clusters for multi period analysis}
\item{lower}{ Lower boundary for the optimal weights (used in \code{\link{optim}})}
\item{upper}{ Upper boundary for the optimal weights (used in \code{\link{optim}})}
\item{...}{ Further arguments to or from other methods.}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
# load examples portfolio
data(ex_ptf)
# set initial weights
ww = c(0.09, rep(0.13,6))
# single period style analysis
Styles(FUND=ex_ptf[,1]
	, IND=ex_ptf[,-1]
	, W=ww
	, lower=NULL
	, upper=NULL)
# multi period style analysis
Multi.Styles(FUND=ex_ptf[,1]
	, IND=ex_ptf[,-1]
	, n_clust=5
	, W=ww
	, lower=NULL
	, upper=NULL)

}




