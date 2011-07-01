\name{Univar}
\alias{univar}
\title{Univariate analysis}
\description{Perform univariate analisys of the dependent variable Y versus each independent
variable X, plotting the results}
\usage{univar(Y, X, stress.period.idx = c(), 
Y.logit = FALSE, Y.logit.adj = 5e-05, 
theme.params = getCurrentTheme(), 
plot = TRUE, overrides = list(...), ...) }
\arguments{ 
\item{Y}{ serie of the dependent variable}
\item{X}{ Matrix containing all independent variables  (one column per variable)}
\item{stress.period.idx}{ vector of positions specifing the stress regime. If provided, the system
will run a modified LS to capture the two regimes}
\item{Y.logit}{ LOGICAL. If TRUE, the dependent variable is transformed using the Logit transform.
The results are retransformed using the inverse Logit. (DEFAULT: FALSE)}
\item{Y.logit.adj}{ Cut-off value. The range of the Y variable is restricted within the interval 
[Y.logit.adj, 1-Y.logit.adj]   (DEFAULT: 0.00005)}
\item{theme.params}{ Theme parameters (DEFAULT: getCurrentTheme())}
\item{plot}{ list of parameters to override the theme. Must match by name the parameters defined by
the theme (DEFAULT: NULL)}
\item{overrides}{ LOGICAL. If TRUE, results are plotted.}
\item{...}{ Further arguments to or from other methods}
}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{ \code{\link{plot.univar}}, \code{\link{print.univar}}  }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
