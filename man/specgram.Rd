\name{specgram}

\alias{specgram}

\title{Spectrogram using short-time Fourier transform}

\description{Computes FFT on each column of X. For Financial series objects (class 'fs'), Close data is extracted.}

\usage{
specgram(X, win.size = max(1, NROW(X)/20), plot = TRUE, ...)
}

\arguments{ 
\item{X}{Matrix of data series (one column per variable).}
\item{win.size}{The size of the window used to compute the FFT}
\item{plot}{LOGICAL. If TRUE, spectrogram is plotted.}
\item{...}{Additional parameters passed to splitWindow, FFT and plot.specgram}
}

\value{
An object of the class 'specgram'. This is an array with dimensions (NFFT, Nwindows, NColX):
\item{NFFT}{The FFT length. It is the next power of 2 greater than the length of each segment/window of X.}
\item{Nwindows}{The number of window segments computed. It depends on the 'by' parameter (default is 1) of the splitWindow function (see details).}
\item{NColX}{The number of columns of X.}

The following attributes are attached to the object:
\item{Fs}{The input Fs parameter to the FFT.}
\item{window}{The window function used to smooth the input data.}
\item{freq}{The frequencies where the FFT was evaluated.}
\item{fpoints}{The array indices where the frequency points relative to 'freq' are stored.}
\item{half}{The input half parameter to the FFT.}
}

\details{
A forward sliding window of length win.size is used to split the input data into segments, then for each segment the FFT of size NFFT = 2^ceiling(log2(win.size)) is computed.
The sliding of the window is controlled by the 'by' parameter of the splitWindow function (default: by = 1). The 'by' parameter should take values between 1 and win.size.
If by = win.size, the input data is split into Nwindows = ceiling(NRowX/win.size) non-overlapping adjacent blocks.
If by = 1 then Nwindows = NRowX - win.size + 1 overlapping segments are computed.
}

\seealso{
\code{\link{splitWindow}}, \code{\link{FFT}}, \code{\link{plot.specgram}}.
}

\examples{
# Load sample financial series data
data(ex_fs)

# 3D spectrogram
specgram(ex_fs, plot3d = TRUE)

# Sampling period
Ts = 0.01
# Generate 10 seconds timeline
t = seq(0, 10, by = Ts)
# Sampling frequency
Fs = 1/Ts
# Linear increasing frequency
f = 2*t
#Chirp signal - Cosine of increasing frequency
chirp = as.matrix(cos(2*pi*f*t))
colnames(chirp) = "Chirp"

# 2D spectrogram
specgram(chirp, Fs = Fs)
# 2D spectrogram with non overlapping windows
specgram(chirp, Fs = Fs, win.size = 128, by = 128)
# 3D spectrogram
specgram(chirp, Fs = Fs, win.size = 128, plot3d = TRUE)
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
