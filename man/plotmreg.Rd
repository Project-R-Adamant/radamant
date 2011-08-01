\name{plotmreg}

\alias{plot.reg}
\alias{plot.mreg}

\title{Plot (Multi)-Regression object}
\description{Plot method for classes 'reg' and 'mreg'.}

\usage{
\method{plot}{mreg}(x, ...) 

\method{plot}{reg}(x
	, mode = c("response", "link")
	, title = ifelse(x$model.type == "lm"
					, "LS Regression"
					, "GLM Regression"
					)
	, theme.params = getCurrentTheme()
	, overrides = list(...)
	, ...
	)
}

\arguments{
\item{x}{Instance of class 'reg'/'mreg'.}
\item{mode}{One of 'response' or 'link'. Controls on which scale results are plotted. See \code{\link{mreg}} for details.}
\item{title}{The plot title}
\item{theme.params}{RAdamant graphics theme.}
\item{overrides}{List of parameters to override the theme.}
\item{...}{Additional arguments passed to cplot.}
}

\value{Void}

\examples{

# Generate some random data
N = 50;
sigma = 1;
X1 = cumsum(rnorm(N));
X2 = cumsum(rnorm(N));

# Define a linear model
Y = 1.5 + X1 + 2*X2 + rnorm(N, sd = sigma);

# Run Multi-Regression
mod = mreg(Y, X = cbind(X1, X2), plot = FALSE);
plot(mod, theme.params = getTheme(2), xlab.srt = 0)

}

\seealso{\code{\link{mreg}}, \code{\link{cplot}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}


