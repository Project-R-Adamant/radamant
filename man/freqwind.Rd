\name{freqwind}

\alias{barthann}
\alias{bartlet}
\alias{blackman}
\alias{gauss}
\alias{lanczos}
\alias{triangle}
\alias{hamming}
\alias{hann}
\alias{cosine}
\alias{kaiser}

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
barthann(N, normalized = TRUE, alpha = 0.38)  
bartlet(N, normalized = TRUE) 
blackman(N, normalized = TRUE, alpha = 0.16) 
cosine(N, normalized = TRUE) 
gauss(N, normalized = TRUE, sigma = 0.5)  
hamming(N, normalized = TRUE) 
hann(N, normalized = TRUE) 
kaiser(N, normalized = TRUE, alpha = 3) 
lanczos(N, normalized = TRUE) 
triangle(N, normalized = TRUE) 
}
\arguments{ 
\item{N}{ Window length}
\item{normalized}{ LOGICAL. If TRUE, window is normalised to have unitary norm}
\item{alpha}{ alpha}
\item{sigma}{ sigma}
}

\value{The window samples}
\author{RAdamant Development Team \email{team@r-adamant.org} }
\seealso{ \code{\link{FFT}} }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
