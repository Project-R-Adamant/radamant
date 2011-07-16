\name{lanczos}

\alias{lanczos}

\title{Lanczos window}

\description{Computes Lanczos window of given length}

\usage{
lanczos(N, normalized = TRUE)
}

\arguments{ 
\item{N}{Window length.}
\item{normalized}{LOGICAL. If TRUE (default), window is normalised to have unitary norm.}
}

\value{An object of the class 'Window'. It is a simple sequence of N samples of the Lanczos window.}

\examples{
\dontrun{
# Generate a Normalised Lanczos window of size 100
x = lanczos(100)
# Plot the window
cplot(x 
	, main = "Lanczos Window" 
	, legend = attr(x, "type") 
	)
# Generate a non-normalised window
y = lanczos(100, FALSE)
# Compare the two
cplot(cbind(x, y) 
	, main = "Lanczos Window" 
	, legend = paste(attr(x, "type"), c("Normalised", "Not Normalised")) 
	, type = c("l", "o") 
	, xlab.srt = 0 
	)
}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
