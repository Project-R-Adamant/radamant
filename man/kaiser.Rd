\name{kaiser}

\alias{kaiser}

\title{Kaiser window}

\description{Computes Kaiser window of given length (Discrete Prolate Spheroidal Sequence approximation).}

\usage{
kaiser(N, normalized = TRUE, alpha = 3)
}

\arguments{ 
\item{N}{Window length.}
\item{normalised}{LOGICAL. If TRUE (default), window is normalised to have unitary norm.}
\item{alpha}{Shape factor (DEFAULT = 3).}
}

\value{An object of the class 'Window'. It is a simple sequence of N samples of the Kaiser window.}

\examples{
# Generate a Kaiser window of size 100
x = kaiser(100, FALSE)
# Plot the window
cplot(x 
	, main = "Kaiser Window" 
	, legend = attr(x, "type") 
	)

# Generate another window with different smoothing factor
y = kaiser(100, normalized = FALSE, alpha = 6)
# Compare the two windows
cplot(cbind(x, y) 
	, main = "Kaiser Window" 
	, legend = paste("Kaiser (alpha = ", c(3, 6), ")", sep = "")
	, type = c("l", "o") 
	, xlab.srt = 0 
	)

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
