\name{Freq_wind}
\alias{barthann}
\alias{bartlet}
\alias{blackman}
\alias{gauss}
\alias{lanczos}
\alias{triangle}
\alias{hamming}
\alias{hann}
\alias{cosine}



\title{Smoothing windows}

\description{
\code{barthann}: Compute Barthann window \cr
\code{bartlet}: Compute Bartlet window \cr
\code{blackman}: Compute Blackman window \cr
\code{gauss}: Compute Gauss window \cr
\code{lanczos}: Compute Lanczos window \cr
\code{triangle}: Compute Triange window \cr
\code{hamming}: Compute Hamming window \cr
\code{hann}: Compute Hann window \cr
\code{cosine}: Compute Cosine window
}

\usage{
barthann(N, normalized = TRUE, alpha = 0.38)  \cr
bartlet(N, normalized = TRUE) \cr
blackman(N, normalized = TRUE, alpha = 0.16) \cr
cosine(N, normalized = TRUE) \cr
gauss(N, normalized = TRUE, sigma = 0.5)  \cr
hamming(N, normalized = TRUE) \cr
hann(N, normalized = TRUE) \cr
kaiser(N, normalized = TRUE, alpha = 3) \cr
lanczos(N, normalized = TRUE) \cr
triangle(N, normalized = TRUE) 
}

\arguments{ 
\item{N}{ Window length}
\item{normalized}{ LOGICAL. If TRUE, window is normalised to have unitary norm}
\item{alpha}{ }
}

\details{}

\value{The window samples}

\references{}

\author{RAdamant Development Team \email{team@r-adamant.org} }

\note{}

\seealso{ \code{\link{FFT}} }

\examples{}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
