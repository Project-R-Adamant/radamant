\name{mtunivar}

\alias{plot.univar}
\alias{print.univar}
\alias{summary.univar}

\title{Methods for univariate analysis}

\description{Print, Plot and Summary methods for class 'univar'}

\usage{
\method{summary}{univar}(object, ...) 

\method{print}{univar}(x, ...) 

\method{plot}{univar}(x, theme.params = getCurrentTheme(), overrides = list(...), ...)
}

\arguments{ 
\item{x, object}{Instance of class 'univar'}
\item{theme.params}{params: Theme parameters (DEFAULT: getCurrentTheme())}
\item{overrides}{List of parameters to override the theme. Must match by name the parameters defined by the theme (DEFAULT: list(...))}
\item{...}{Alternative way to quickly override theme parameters.}
}

\examples{
# Load sample time series data
data(ex_ptf)
# Define the dependent variable
Y = ex_ptf[, 1, drop = FALSE];
# Define the independent variables
X = ex_ptf[, -1];
# Define x-axis labels
time.labels = paste("t[", 1:length(Y), "]", sep = "")
# Univar analysis
res = univar(Y, X, plot = FALSE);
plot(res
    , xlabels = parse(text = time.labels)
    # Remove x-labels rotation
    , xlab.srt = 0
    # Set more space between x-labels and the x-axis line (10% of diff(par("usr")[3:4]))
    , xlab.offset = 0.1
    # Set more space between x-title and the x-axis line  (20% of diff(par("usr")[3:4]))
    , xtitle.offset = 0.2
    # Only 4 tickmarks on the y-axis
    , y.ticks = 4
    )
}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{univar}}, \code{\link{cplot}}}


