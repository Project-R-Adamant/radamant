\name{resvecar}

\alias{resid.VecAr}
\alias{residuals.VecAr}
\alias{resid.VecAr}
\alias{residuals.VecAr}

\title{Extract Model Residuals from Vector AutoRegressive object}
\description{Generic method for extracting model residuals from object of class 'VecAr'.}

\usage{
\method{residuals}{VecAr}(object, na.rm = FALSE, ...)
}

\arguments{
\item{object}{Instance of class 'VecAr'.} 
\item{na.rm}{Logical. If TRUE, NA records are removed.}
\item{...}{Further arguments to or from other methods.}
}

\value{A matrix containing all model residuals, one column for each model.}


\examples{
# Collect series data
X = cbind(BJsales, BJsales.lead);

# Generate simple VAR(2) model
mod = VecAr(X, ar.lags = 1:2)

# Extract residuals (note NA due to the lagged data)
residuals(mod)
residuals(mod, na.rm = TRUE)
}

\seealso{\code{\link{VecAr}}, \code{\link{residuals.mreg}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}


