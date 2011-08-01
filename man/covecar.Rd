\name{covecar}

\alias{coef.VecAr}

\title{Extract Model Coefficients from Vector AutoRegressive object}
\description{Generic method for extracting model coefficients matrix from object of class 'VecAr'.}

\usage{
\method{coef}{VecAr}(object, ...)
}

\arguments{
\item{object}{Instance of class 'VecAr'.} 
\item{...}{Further arguments to or from other methods.}
}

\value{A matrix containing all model coefficients, one column for each variable in the model.}

\examples{
# Collect series data
X = cbind(BJsales, BJsales.lead);

# Generate simple VAR(2) model
mod = VecAr(X, ar.lags = 1:2)

# Extract coefficients
coef(mod)
}

\seealso{\code{\link{VecAr}}, \code{\link{coef.mreg}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}


