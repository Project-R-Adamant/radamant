\name{plotfft}

\alias{plot.FFT}

\title{Customised Fast Fourier Transform - Plotting}

\description{Plot function for class 'FFT'. Plots Modulus and Phase for each column of the FFT object x}

\usage{

\method{plot}{FFT}(x 
	, theme.params = getCurrentTheme() 
	, overrides = list(...) 
	, shaded = TRUE 
	, show.periodicity = FALSE 
	, show.legend = FALSE 
	, zoom = 100 
	, semilog = FALSE 
	, new.device = FALSE 
	, ... 
	) 

}

\arguments{ 
\item{x}{Instance of class 'FFT'.}
\item{theme.params}{theme parameters list (DEFAULT: getCurrentTheme()).}
\item{overrides}{List of parameters to override the theme. Only parameters that match those defined by the theme are overridden (DEFAULT: list(...)).}
\item{shaded}{LOGICAL. If TRUE, the modulus of x is shaded.}
\item{show.periodicity}{LOGICAL. If TRUE, Periods (1/frequencies) are showed instead of frequencies on the x-axis (DEFAULT = FALSE).}
\item{show.legend}{LOGICAL. If TRUE, legend is added to the plot (DEFAULT = FALSE)}
\item{...}{Additional parameters passed to the cplot function. Also used to quickly specify theme overrides.}
}

\value{Void}

\seealso{
\code{\link{cplot}}.
}

\examples{
# Load sample financial series data
data(ex_fs)

# Frequency Analysis
Xf = FFT(ex_fs, plot = FALSE)

# Plot full spectrum
plot(Xf)

# Plot falf spectrum (right side) and use blackman windowing, remove area shading
plot(Xf, half = TRUE, window = blackman, shaded = FALSE)

# Show periodicity instead of frequency, and use hamming window
plot(Xf, half = TRUE, window = hamming, show.periodicity = TRUE)

# Use kaiser window, zoom in to show only 10% of the half frequency spectrum, use semilog axis (decibel)
plot(Xf, half = TRUE, window = kaiser, show.periodicity = TRUE, zoom = 10, semilog = TRUE)

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
