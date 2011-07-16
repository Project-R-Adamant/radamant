\name{blackman}

\alias{blackman}

\title{Blackman window}

\description{Computes Blackman window of given length}

\usage{
blackman(N, normalized = TRUE, alpha = 0.16)
}

\arguments{ 
\item{N}{Window length.}
\item{normalised}{LOGICAL. If TRUE (default), window is normalised to have unitary norm.}
\item{alpha}{Shape factor (DEFAULT = 0.16). Determines the smoothing of the window's sidelobes.}
}

\value{An object of the class 'Window'. It is a simple sequence of N samples of the Blackman window.}

\examples{
# Generate a Blackman window of size 100
x = blackman(100, FALSE)
# Plot the window
cplot(x 
	, main = "Blackman Window" 
	, legend = attr(x, "type") 
	)

# Generate another window with lower smoothing factor
y = blackman(100, normalized = FALSE, alpha = 0.4)
# Compare the two windows
cplot(cbind(x, y) 
	, main = "Blackman Window" 
	, legend = paste("Blackman (alpha = ", c(0.16, 0.4), ")", sep = "")
	, type = c("l", "o") 
	, xlab.srt = 0 
	)

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
