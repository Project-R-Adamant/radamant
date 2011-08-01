\name{univar}
\alias{univar}
\title{Univariate analysis}
\description{Perform univariate analisys of the dependent variable Y versus each independent variable X, plotting the results}
\usage{
univar(Y
    , X
    , stress.period.idx = c()
    , Y.logit = FALSE
    , Y.logit.adj = 0.00005
    , theme.params = getCurrentTheme()
    , plot = TRUE
    , overrides = list(...)
    , ...
    )
}
\arguments{ 
\item{Y}{The dependent variable. This must be a one column matrix.}
\item{X}{Matrix containing all independent variables  (one column per variable)}
\item{stress.period.idx}{Vector of positions specifing the stress regime. 
If provided, the system will run a modified LS to capture the two regimes}
\item{Y.logit}{Logical. If TRUE, the dependent variable is transformed using the Logit transform. 
Results are then transformed back using the inverse Logit. (Default: FALSE)}
\item{Y.logit.adj}{Cut-off value. The range of the Y variable is restricted within the interval [Y.logit.adj, 1-Y.logit.adj] (Default: 0.00005)}
\item{theme.params}{Theme parameters (Default: getCurrentTheme())}
\item{plot}{Logical. If TRUE, results are plotted.}
\item{overrides}{List of parameters to override the theme. Must match by name the parameters defined by the theme (Default: list(...)).}
\item{...}{Alternative way to quickly override theme parameters.}
}

\value{
An object of class 'univar'. This is a list with the following components:
\item{Y.logit}{The input Y.logit parameter.}
\item{stress.idx}{The input stress.period.idx parameter.}
\item{model}{A list of NCOL(X) entries. Each entry is a linear model object (of class 'lm'): regression Y on the corresponding column of X.}
\item{summary}{A summary data frame with columns [regressor, formula, eq, sigma.squared, adj.r.squared, pvalue].}
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
# Univar Analysis
univar(Y, X
    , xlabels = parse(text = time.labels)
    # Remove x-labels rotation
    , xlab.srt = 0
    # Set more space between x-labels and the x-axis line 
        # (10% of diff(par("usr")[3:4]))
    , xlab.offset = 0.1
    # Set more space between x-title and the x-axis line  
        # (20% of diff(par("usr")[3:4]))
    , xtitle.offset = 0.2
    # Only 4 tickmarks on the y-axis
    , y.ticks = 4
    )
}


\author{RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{plot.univar}}, \code{\link{print.univar}} }
