\name{covesvar}

\alias{estVar.VecAr}
\alias{vcov.VecAr}

\title{Compute residual and coefficients covariance matrix from Vector AutoRegressive object}

\description{Generic method for computing residual and coefficients covariance matrix from object of class 'VecAr'.}

\usage{
\method{estVar}{VecAr}(object, ...)
\method{vcov}{VecAr}(object, ...)
}

\arguments{
\item{object}{Instance of class 'VecAr'.} 
\item{...}{Further arguments to or from other methods.}
}

\value{A matrix with calculated residual / coefficients covariance}

\examples{
# Collect series data
X = cbind(BJsales, BJsales.lead);

# Generate simple VAR(2) model
mod = VecAr(X, ar.lags = 1:2)

# Extract residual covariance matrix
estVar.VecAr(mod)
# Extract coefficients covariance matrix
vcov.VecAr(mod)

}

\seealso{\code{\link{VecAr}}, \code{\link{residuals.VecAr}}, \code{\link{coef.VecAr}} }

\author{RAdamant Development Team \email{team@r-adamant.org}}


