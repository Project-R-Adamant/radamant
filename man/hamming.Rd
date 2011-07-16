\name{hamming}

\alias{hamming}

\title{Hamming window}

\description{Computes Hamming window of given length}

\usage{
hamming(N, normalized = TRUE)
}

\arguments{ 
\item{N}{Window length.}
\item{normalised}{LOGICAL. If TRUE (default), window is normalised to have unitary norm.}
}

\value{An object of the class 'Window'. It is a simple sequence of N samples of the Hamming window.}

\examples{
# Generate a Normalised Hamming window of size 100
x = hamming(100)
# Plot the window
cplot(x 
	, main = "Hamming Window" 
	, legend = attr(x, "type") 
	)
# Generate a non-normalised window
y = hamming(100, FALSE)
# Compare the two
cplot(cbind(x, y) 
	, main = "Hamming Window" 
	, legend = paste(attr(x, "type"), c("Normalised", "Not Normalised")) 
	, type = c("l", "o") 
	, xlab.srt = 0 
	)

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
