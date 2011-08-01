\name{sensplot}

\alias{plot.sensAnalysis}

\title{Plot Sensitivity Analysis}

\description{Plot method for class 'sensAnalysis'.}

\usage{
\method{plot}{sensAnalysis}(x
	, main = NULL
	, xlabels = rownames(x$coeffs)
	, xtitle = ""
	, theme.params = getCurrentTheme()
	, ...
	)
}

\arguments{
\item{x}{A Sensitivity Analysis object. Instance of class 'sensAnalysis'.}
\item{main}{Main plot title}
\item{xlabels}{Labels for the x-axis}
\item{xtitle}{Title for the x-axis}
\item{theme.params}{RAdamant graphics theme.}
\item{...}{Further arguments passed to the cplot function.}
}

\value{Void}

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
res = sensAnalysis(mod
		# Sliding Window with 20 samples
		, mode = "SW"
		, win.size = 20
		# Shift by 5 points backward at each step
		, direction = "backward"
		, by = 5
		);

# Plot results for the first model
plot(res[[1]]
	# Use another theme
	, theme.params = getTheme(2)
	# Override theme - show all labels on the x-axis
	, x.ticks = "ALL"
	)

}

\seealso{\code{\link{sensAnalysis}}, \code{\link{mreg}}, \code{\link{plot.sensAnalysis}}, \code{\link{cplot}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}
