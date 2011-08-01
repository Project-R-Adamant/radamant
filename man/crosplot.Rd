\name{crosplot}

\alias{cross.plot}
\title{Y Vs X Cross Plot}

\description{Plot the input dependent variable Y versus each input independent variable X}

\usage{
cross.plot(Y
    , X
    , theme.params = getCurrentTheme()
    , xlabels = NULL
    , two.axis = TRUE
    , shaded.first = FALSE
    , overrides = list(...)
    , ...
    )
}

\arguments{ 
\item{Y}{Dependent variable.}
\item{X}{Matrix containing all independent variables (one column per variable).}
\item{theme.params}{Theme parameters (Default: getCurrentTheme()).}
\item{xlabels}{Vector of labels associated to the rows of X (i.e. Time labels)(Default: NULL)}
\item{two.axis}{Logical. If TRUE, series are plotted on two axis (two scales).}
\item{shaded.first}{Logical. If TRUE, the variable Y is shaded.}
\item{overrides}{List of parameters to override the theme. Must match by name the parameters defined by the theme (Default: list(...))}
\item{...}{Alternative way to quickly override the theme.}
}

\value{Void}

\examples{
# Load sample time series data
data(ex_ptf)
# Define the dependent variable
Y = ex_ptf[, 1, drop = FALSE];
# Define the independent variables
X = ex_ptf[, -1];
# Define x-axis labels
time.labels = paste("t[", 1:length(Y), "]", sep = "")
# Cross plot
cross.plot(Y, X
        , xlabels = parse(text = time.labels)
        , overrides = list(xlab.srt = 0)
        )
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
