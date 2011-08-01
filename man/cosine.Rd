\name{cosine}

\alias{cosine}

\title{Cosine window}

\description{Computes Cosine window of given length}

\usage{
cosine(N, normalized = TRUE)
}

\arguments{ 
\item{N}{Window length.}
\item{normalized}{Logical. If TRUE (default), window is normalised to have unitary norm.}
}

\value{An object of the class 'Window'. It is a simple sequence of N samples of the Cosine window.}

\examples{
# Generate a Normalised Cosine window of size 100
x = cosine(100)
# Plot the window
cplot(x 
    , main = "Cosine Window" 
    , legend = attr(x, "type") 
    )
# Generate a non-normalised window
y = cosine(100, FALSE)
# Compare the two
cplot(cbind(x, y) 
    , main = "Cosine Window" 
    , legend = paste(attr(x, "type"), c("Normalised", "Not Normalised")) 
    , type = c("l", "o") 
    , xlab.srt = 0 
    )

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
