\name{bartlet}

\alias{bartlet}

\title{Bartlet window}

\description{Computes Bartlet window of given length}

\usage{
bartlet(N, normalized = TRUE)
}

\arguments{ 
\item{N}{Window length.}
\item{normalized}{LOGICAL. If TRUE (default), window is normalised to have unitary norm.}
}

\value{An object of the class 'Window'. It is a simple sequence of N samples of the Bartlet window.}

\examples{
# Generate a Normalised Bartlet window of size 100
x = bartlet(100)
# Plot the window
cplot(x 
	, main = "Bartlet Window" 
	, legend = attr(x, "type") 
	)
# Generate a non-normalised window
y = bartlet(100, FALSE)
# Compare the two
cplot(cbind(x, y) 
	, main = "Bartlet Window" 
	, legend = paste(attr(x, "type"), c("Normalised", "Not Normalised")) 
	, type = c("l", "o") 
	, xlab.srt = 0 
	)

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
