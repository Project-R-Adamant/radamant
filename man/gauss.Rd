\name{gauss}
\alias{gauss}
\title{Gauss window}
\description{Computes Gauss window of given length}
\usage{
gauss(N, normalized = TRUE, sigma = 0.5)
}
\arguments{ 
\item{N}{Window length.}
\item{normalized}{LOGICAL. If TRUE (default), window is normalised to have unitary norm.}
\item{sigma}{Standard Deviation - Expansion factor. sigma <= 0.5.}
}
\value{An object of the class 'Window'. It is a simple sequence of N samples of the Gauss window.}
\examples{
# Generate a Normalised Gauss window of size 100
x = gauss(100)
# Plot the window
cplot(x 
	, main = "Gauss Window" 
	, legend = attr(x, "type") 
	)
	
# Generate a non-normalised window
y = gauss(100, FALSE)
# Compare the two
cplot(cbind(x, y) 
	, main = "Gauss Window" 
	, legend = paste(attr(x, "type"), c("Normalised", "Not Normalised")) 
	, type = c("l", "o") 
	, xlab.srt = 0 
	)

# Generate another window with smaller expansion factor
z = gauss(100, normalized = FALSE, sigma = 0.1)
# Compare the two expansion factors
cplot(cbind(y, z) 
	, main = "Gauss Window" 
	, legend = paste("Gauss (sigma = ", c(0.5, 0.1), ")")
	, type = c("l", "o") 
	, xlab.srt = 0 
	)
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
