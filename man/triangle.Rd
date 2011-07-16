\name{triangle}

\alias{triangle}

\title{Triangle window}

\description{Computes Triangle window of given length}

\usage{
triangle(N, normalized = TRUE)
}

\arguments{ 
\item{N}{Window length.}
\item{normalised}{LOGICAL. If TRUE (default), window is normalised to have unitary norm.}
}

\value{An object of the class 'Window'. It is a simple sequence of N samples of the Triangle window.}

\examples{
# Generate a Normalised Triangle window of size 100
x = triangle(100)
# Plot the window
cplot(x 
	, main = "Triangle Window" 
	, legend = attr(x, "type") 
	)
# Generate a non-normalised window
y = triangle(100, FALSE)
# Compare the two
cplot(cbind(x, y) 
	, main = "Triangle Window" 
	, legend = paste(attr(x, "type"), c("Normalised", "Not Normalised")) 
	, type = c("l", "o") 
	, xlab.srt = 0 
	)

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
