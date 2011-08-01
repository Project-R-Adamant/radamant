\name{plotspec}

\alias{plot.specgram}

\title{Spectrogram Plotting}

\description{Plot method for class 'specgram'.}

\usage{
\method{plot}{specgram}(x 
    , show.periodicity = FALSE 
    , theme.params = getCurrentTheme() 
    , xtitle = "Time" 
    , ytitle = ifelse(show.periodicity
                     , "Periodicity"
                     , "Frequency"
                     ) 
    , plot3d = FALSE 
    , overrides = list(...) 
    , ...
    )
}

\arguments{ 
\item{x}{Instance of class 'specgram'}
\item{show.periodicity}{Logical. If TRUE, Periods (1/frequencies) are showed instead of frequencies on the x-axis (Default: FALSE)}
\item{theme.params}{RAdamant graphics theme. (Default: getCurrentTheme())}
\item{xtitle}{Title for the x-axis (Default: "Time")}
\item{ytitle}{Title for the y-axis (Default: "Frequency" or "Periodicity" depending on the value of show.periodicity)}
\item{plot3d}{Logical. If TRUE, 3D spectrogram is plotted.}
\item{overrides}{List of parameters to override the theme. Only parameters that match those defined by the theme are overridden (Default: list(...))}
\item{...}{Used to quickly specify theme overrides.}
}

\value{
Void
}


\seealso{
\code{\link{specgram}}.
}

\examples{
# Load sample financial series data
data(ex_fs)

# 3D spectrogram
spec = specgram(ex_fs, plot = FALSE)
# Plotting
plot(spec, plot3d = TRUE)
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
