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
\item{theme.params}{RAdamant graphics theme (Default: getCurrentTheme()).}
\item{overrides}{List of parameters to override the theme. Only parameters that match those defined by the theme are overridden (Default: list(...)).}
\item{shaded}{Logical. If TRUE, the modulus of x is shaded.}
\item{show.periodicity}{Logical. If TRUE, Periods (1/frequencies) are showed instead of frequencies on the x-axis (Default: FALSE).}
\item{show.legend}{Logical. If TRUE, legend is added to the plot (Default: FALSE).}
\item{zoom}{ Zoom}
\item{semilog}{Logical. If TRUE, the modulus of the FFT is shown on a dB scale.}
\item{new.device}{Logical. If TRUE, a new plotting device is opened.}
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

# Plot half spectrum
plot(Xf
    # Plot half spectrum (right side)
    , half = TRUE
    # Use Blackman window
    , window = blackman
    # Remove area shading
    , shaded = FALSE)

# Show periodicity instead of frequency, and use hamming window
plot(Xf
    , half = TRUE
    , window = hamming
    , show.periodicity = TRUE)

# Plot with other options 
plot(Xf
    # Half spectrum
    , half = TRUE
    # Kaiser Window
    , window = kaiser
    , # Show cycles instead of frequencies
    , show.periodicity = TRUE
    # Zoom in to show only 10% of the half spectrum
    , zoom = 10
    # Use semilog axis (decibel)
    , semilog = TRUE
    )

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
