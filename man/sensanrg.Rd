\name{sensanrg}

\alias{sensAnalysis.reg}
\alias{sensAnalysis.mreg}

\title{Sensitivity Analysis for Multi-Regression Models}

\description{Sensitivity analysis method for classes 'reg' and 'mreg'.}

\usage{
\method{sensAnalysis}{reg}(X, ...)

\method{sensAnalysis}{mreg}(X, ...)
}

\arguments{
\item{X}{A regression model. Instance of class 'reg', 'mreg'.}
\item{...}{Further arguments passed to the default method.}
}

\value{
An instance of class 'sensAnalysis' if X has class 'reg', or a list of length(X) objects of class 'sensAnalysis' if X has class 'mreg'.
}

\examples{
# Generate some random data
N = 50;
sigma = 0.1;
X1 = cumsum(rnorm(N));
X2 = rnorm(N);
X3 = cumsum(rnorm(N));
X4 = rnorm(N);

# Define a linear model
Y1 = 1.5 + X1 + 2*X3 + rnorm(N, sd = sigma);
# Define a logit model
Y2 = inv.logit(-2.2 + 0.3*X2 - 0.2*X4 + rnorm(N, sd = sigma));

# Run Multi-Regression
mod = mreg(Y = cbind(Y1, Y2)
		, X = cbind(X1, X2, X3, X4)
		# Stepwise regression
		, type = "stepwise"
		# lm on Y1 and glm on Y2
		, mode = c("lm", "glm")
		# Set the family. 
		# It is recycled but family is only used for glm
		, family = "binomial"
		# Constrain the maximum number of variables 
		# that can enter the regression
		, max.vars = c(3, 2)
		# Use another theme
		, theme.params = getTheme(2)
		);

# Perform Sensitivity Analysis, Backward Sliding Window
sensAnalysis(mod
	# Sliding Window with 20 samples
	, mode = "SW"
	, win.size = 20
	# Shift by 5 points backward at each step
	, direction = "backward"
	, by = 5
	# Plot results
	, plot = TRUE
	# Override theme - show all labels on the x-axis
	, x.ticks = "ALL"
	)

}

\seealso{\code{\link{sensAnalysis}}, \code{\link{mreg}}, \code{\link{plot.sensAnalysis}}, \code{\link{cplot}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}


