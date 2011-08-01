\name{sensan}

\alias{sensAnalysis.default}
\alias{sensAnalysis}

\title{Sensitivity Analysis}

\description{Generic method for parameter sensitivity analysis on regression models.}

\usage{
sensAnalysis(X, ...)

\method{sensAnalysis}{default}(X, win.size = length(coef(X)), plot = FALSE, ...)
}

\arguments{
\item{X}{A regression model. Instance of class 'lm', 'glm'.}
\item{win.size}{The initial window size for the analysis. See \code{\link{splitWindow}} for details.}
\item{plot}{Logical. If TRUE, results are plotted.}
\item{...}{Further arguments passed to splitWindow and cplot.}
}

\value{
An object of class 'sensAnalysis'. This is a list with the following elements:
\item{coeffs}{Matrix of regression coefficients estimated on each portion of data delimited by the indexes computed by splitWindow.}
\item{weights}{Matrix of regression weights as computed by get.lm.weights.}
\item{pvalues}{Matrix of p-values of the regression coefficiens.}
}

\examples{
# Generate some random data
N = 50;
sigma = 1;
X1 = cumsum(rnorm(N));
X2 = cumsum(rnorm(N));

# Define a linear model
Y = 1.5 + X1 + 2*X2 + rnorm(N, sd = sigma);

# Run Regression
mod = lm(Y ~ X1 + X2);

# Perform Sensitivity Analysis, Forward Extended Window (Default)
sensAnalysis(mod
    # Starting with 10 samples
    , win.size = 10
    # Increment by 5 points at each step
    , by = 5
    )

}

\seealso{\code{\link{splitWindow}}, \code{\link{get.lm.weights}}, \code{\link{plot.sensAnalysis}}, \code{\link{cplot}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}
