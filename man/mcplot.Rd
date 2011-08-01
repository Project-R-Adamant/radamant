\name{mcplot}
\alias{mcplot}
\title{Muliple correlation plot}
\description{Multiple correlation plot}

\usage{
mcplot(X
    , hist.nclass = 10
    , theme.params = getCurrentTheme()
    , coLin = TRUE
    , main = ifelse(coLin, "Co-Linearity Analysis"
                         , "Multi-Correlation Analysis")
    , new.device = FALSE
    , ...
    ) 
}

\arguments{
\item{X}{Matrix of data series (one column per variable).}
\item{hist.nclass}{Number of bins used for computing histogram plot (Default: 10).}
\item{theme.params}{RAdamant graphics theme.}
\item{coLin}{Logical. If TRUE, Co-Linearity analysis is performed, otherwise Correlation analysis is assumed. See details.}
\item{main}{The plot title}
\item{new.device}{Logical. If TRUE, a new devide is opened.}
\item{...}{ Further arguments passed to chist.}
}

\value{Void}

\details{
The parameter 'coLin' controls how correlation coefficients are displayed:
\itemize{
\item{coLin = TRUE: the higher the correlation (in absolute terms) the more the corresponding columns are collinear. \cr
The correlation coefficient is displayed with variable colors ranging from green (abs(rho) = 0) to red (abs(rho) = 1).}

\item{coLin = FALSE: Colors are switched ranging from red (abs(rho) = 0) to green (abs(rho) = 1).}
}
}

\examples{
# Load sample time series data
data(ex_ptf);

# Compute Multi Collinearity Analysis 
# (High correlation (abs(rho)) in red)
mcplot(ex_ptf[, c(2:5)]);

# Compute Multi Correlation Analysis 
# (High correlation (abs(rho)) in green)
mcplot(ex_ptf[, c(2:5)]
    # Increase number of histogram bins
    , hist.nclass = 30
    # Specify correlation type analysis
    , coLin = FALSE
    # Use Normal distribution fitting for the histograms
    , density = "normal"
    );
}

\seealso{\code{\link{chist}}}

\author{RAdamant Development Team \email{team@r-adamant.org}}
