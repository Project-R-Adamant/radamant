\name{fft}

\alias{FFT}
\alias{FFT.default}

\title{Customised Fast Fourier Transform}

\description{Computes FFT on each column of X. For Financial series objects (class 'fs'), Close data is extracted.}

\usage{
FFT(x, ...)

\method{FFT}{default}(x
    , Fs = 1
    , half = FALSE
    , window = NULL
    , plot = TRUE
    , optimised = TRUE
    , ...
    ) 

}

\arguments{ 
\item{x}{Matrix of data series (one column per variable).}
\item{Fs}{Sampling frequency (Default: 1).}
\item{half}{Logical. If TRUE, half spectrum indices are computed.}
\item{window}{Function or character name of the window used to smooth the data (Default: NULL. Results in rectangular window).}
\item{plot}{Logical. If TRUE, frequency spectrum is plotted.}
\item{optimised}{Logical. If TRUE, the number of FFT evaluation points is the next integer (power of 2) that allows the fast computation}
\item{...}{Additional parameters passed to the plot (in the default implementation)}
}

\value{An object of the class 'FFT'. It is a complex matrix (same number of columns as x) of frequency data.
The following attributes are attached to the object:
\item{Fs}{The input Fs parameter}
\item{window}{The window function used to smooth the input data}
\item{freq}{The frequencies where the FFT was evaluated}
\item{fpoints}{The array indices where the frequency points relative to 'freq' are stored}
\item{half}{The input half parameter.}
}

\examples{
# Load sample financial series data
data(ex_fs)

# Frequency Analysis - Full spectrum
FFT(ex_fs)

# Frequency Analysis - Half spectrum (right side) and use blackman windowing, remove area shading
FFT(ex_fs
    , half = TRUE
    , window = blackman
    , shaded = FALSE)

# Show periodicity instead of frequency, and use hamming window
FFT(ex_fs
    , half = TRUE
    , window = hamming
    , show.periodicity = TRUE)

# Use kaiser window, zoom in to show only 10% of the half frequency spectrum, use semilog axis (decibel)
FFT(ex_fs
    , half = TRUE
    , window = kaiser
    , show.periodicity = TRUE
    , zoom = 10
    , semilog = TRUE)

# Multiple FFT on matrix input. 
# Use Bartlet-Hann window, zoom in to show only 20% of the full frequency spectrum, use semilog axis (decibel)
FFT(ex_fs[,]
    , window = barthann
    , zoom = 20
    , semilog = TRUE
    , shaded = FALSE)

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
