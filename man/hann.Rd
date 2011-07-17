\name{hann}

\alias{hann}

\title{Hann window}

\description{Computes Hann window of given length}

\usage{
hann(N, normalized = TRUE)
}

\arguments{ 
\item{N}{Window length.}
\item{normalized}{LOGICAL. If TRUE (default), window is normalised to have unitary norm.}
}

\value{An object of the class 'Window'. It is a simple sequence of N samples of the Hann window.}

\examples{
# Generate a Normalised Hann window of size 100
x = hann(100)
# Plot the window
cplot(x 
	, main = "Hann Window" 
	, legend = attr(x, "type") 
	)
# Generate a non-normalised window
y = hann(100, FALSE)
# Compare the two
cplot(cbind(x, y) 
	, main = "Hann Window" 
	, legend = paste(attr(x, "type"), c("Normalised", "Not Normalised")) 
	, type = c("l", "o") 
	, xlab.srt = 0 
	)

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
