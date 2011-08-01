\name{colinred}
\alias{colin.reduce}
\title{Multi Co-Linearity reduction}
\description{
Performs multicollinearity reduction.
Cross Co-Linearity analysis between the columns of Y and X is computed, then for each column Yi, 
a reduced set of the columns of X is computed by removing those columns that are too correlated (one for each co-linear pair). \cr
In the removal process, those columns of X that are most correlated to Yi are kept.
}

\usage{
colin.reduce(Y, X, max.iter = 100, trsh = 0.85)
}

\arguments{
\item{Y}{Matrix of data series - Dependent variables (one column per variable).}
\item{X}{Matrix of data series - Independent variables (one column per variable).}
\item{max.iter}{Max number of iterations allowed.}
\item{trsh}{Threshold over which two columns are considered too correlated (Default: 0.8).}
}

\value{
A list of Ny elements (Ny = number of columns of Y):
\item{i-th element}{Matrix containing a subset of the columns of X. This is obtained by removing collinear entries. \cr
This element of the list is named after the corresponding i-th column of Y (or a default is given if Yi has no name).}
}

\examples{
# Load sample time series data
data(ex_ptf);

# Select dependent variable
Y = ex_ptf[, 1, drop = FALSE];
# Select indipendent variables
X = ex_ptf[, -1, drop = FALSE];
# Print column names
colnames(X)

# Collinearity Reduction
X.red = colin.reduce(Y, X, trsh = 0.8);
# Print remaining column names
colnames(X.red[[1]])
}

\seealso{\code{\link{colin.pairs}}, \code{\link{cross.colin}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}
