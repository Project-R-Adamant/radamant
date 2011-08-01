\name{colinprs}
\alias{colin.pairs}
\title{Co-Linearity analysis}
\description{This function performs a Co-Linearity analysis between the columns of X.\cr
Correlation factors between columns are computed, and pairs of columns with a correlation factor higher than a specified threshold are returned.
}

\usage{colin.pairs(X, trsh = 0.8)}

\arguments{
\item{X}{Matrix of data series (one column per variable).}
\item{trsh}{Threshold over which two columns are considered too correlated (Default: 0.8).}
}

\value{
A list of with the following elements:
\item{CoLinMat}{Lower Triangular correlation matrix (Correlations between the columns of X).}
\item{CoLinPairs}{Data frame of columns [VAR1, VAR2, Rho] containing the pairs of columns with a correlation factor higher than the given threshold, sorted in descending order.}
}

\examples{
# Load sample time series data
data(ex_ptf);

# Compute correlation matrix and column pairs 
	# with correlation higher than 0.8
colin.pairs(ex_ptf);
}

\author{RAdamant Development Team \email{team@r-adamant.org}}