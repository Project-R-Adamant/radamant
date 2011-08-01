\name{crscolin}
\alias{cross.colin}
\title{Cross Co-Linearity Analysis}
\description{Perform a cross Co-Linearity analysis between the columns of Y and X: \cr
Correlation factors between each column Yi and all columns of X are calculated for different time lags. \cr
Pairs of columns of X with a correlation factor higher than a specified threshold are also returned.
}

\usage{cross.colin(Y, X, max.lag = 8, trsh = 0.8)}
\arguments{
\item{Y}{Matrix of data series - Dependent variables (one column per variable)}
\item{X}{Matrix of data series - Independent variables (one column per variable)}
\item{max.lag}{Max lag for which cross correlation is computed}
\item{trsh}{Threshold over which two columns are considered too correlated (Default: 0.8)}
}

\value{
A list of Ny + 2 elements (Ny = number of columns of Y):
\item{First Ny elements}{Lagged correlation matrix (Nx by max.lag+1) between Yi and X. Named as the column names of Y (or default is given if null).}
\item{CoLinMat}{Lower Triangular correlation matrix (Correlations between the columns of X)}
\item{CoLinPairs}{Data frame of columns [VAR1, VAR2, Rho] containing the pairs of columns with a correlation factor higher than the given threshold, sorted by Rho in descending order.}
}

\examples{
# Load sample time series data
data(ex_ptf);

# Select dependent variable
Y = ex_ptf[, 1, drop = FALSE];
# Select indipendent variables
X = ex_ptf[, -1, drop = FALSE];

# Cross Co-Linearity analysis
cross.colin(Y, X, max.lag = 4, trsh = 0.8);
}

\seealso{\code{\link{colin.pairs}}}


\author{RAdamant Development Team \email{team@r-adamant.org}}
