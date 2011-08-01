\name{fitvecar}

\alias{fitted.VecAr}

\title{Extract Model Fitted Values from Vector AutoRegressive object}
\description{Generic method for extracting model fitted values from object of class 'VecAr'.}

\usage{
\method{fitted}{VecAr}(object, ...)
}

\arguments{
\item{object}{Instance of class 'VecAr'.} 
\item{...}{Further arguments to or from other methods.}
}

\examples{
# Collect series data
X = cbind(BJsales, BJsales.lead);

# Generate simple VAR(2) model
mod = VecAr(X, ar.lags = 1:2)

# Extract fitted values
fitted(mod)
}

\seealso{\code{\link{VecAr}}, \code{\link{predict.mreg}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}


