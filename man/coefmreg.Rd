\name{coefmreg}

\alias{coef.reg}
\alias{coef.mreg}

\title{Extract Model Coefficients for (Multi)-Regression object}

\description{Generic method for extracting model coefficoents from object of classes 'reg' and 'mreg'.}

\usage{
\method{coef}{reg}(object, ...)

\method{coef}{mreg}(object, ...)
}

\arguments{ 
\item{object}{Instance of class 'reg'/'mreg'.}
\item{...}{Further arguments to or from other methods.}
}

\value{One of the following:
\itemize{
    \item{class 'mreg': A matrix containing all model coefficients, one column for each model.}
    \item{class 'reg': A matrix containing the model specific coefficients.}
}
}

\seealso{\code{\link{mreg}}.}

\examples{
# Generate some random data
N = 50;
sigma = 1;
X1 = cumsum(rnorm(N));
X2 = cumsum(rnorm(N));

# Define a linear model
Y1 = 1.5 + X1 + 2*X2 + rnorm(N, sd = sigma);
Y2 = -2 + 1.2*X1 -X2 + rnorm(N, sd = sigma);

# Run Multi-Regression
mod = mreg(Y = cbind(Y1, Y2), X = cbind(X1, X2), plot = FALSE);
# Extract all coefficients
coef(mod)
# Extract coefficients from the first model
coef(mod[[1]])
}


\author{RAdamant Development Team \email{team@r-adamant.org}}


