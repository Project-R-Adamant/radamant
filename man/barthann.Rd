\name{barthann}

\alias{barthann}

\title{Bartlet-Hann window}

\description{Computes Bartlet-Hann window of given length}

\usage{
barthann(N, normalized = TRUE, alpha = 0.38)
}

\arguments{ 
\item{N}{Window length.}
\item{normalised}{LOGICAL. If TRUE (default), window is normalised to have unitary norm.}
\item{alpha}{Shape factor (DEFAULT = 0.38).}
}

\value{An object of the class 'Window'. It is a simple sequence of N samples of the Bartlet-Hann window.}

\examples{
# Generate a Bartlet-Hann window of size 100
x = barthann(100, FALSE)
# Plot the window
cplot(x 
	, main = "Bartlet-Hann Window" 
	, legend = attr(x, "type") 
	)

# Generate another window with different smoothing factor
y = barthann(100, normalized = FALSE, alpha = 0.5)
# Compare the two windows
cplot(cbind(x, y) 
	, main = "Bartlet-Hann Window" 
	, legend = paste("Bartlet-Hann (alpha = ", c(0.38, 0.5), ")", sep = "")
	, type = c("l", "o") 
	, xlab.srt = 0 
	)

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
