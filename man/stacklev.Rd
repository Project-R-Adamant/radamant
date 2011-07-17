\name{stacklev}
\alias{CallStackLevels}
\title{Retrieve the number of calls in the stack.}
\description{
Retrieve the number of calls in the stack. To be called from inside a function.
}

\usage{
CallStackLevels()
}

\value{
The number of calls in the stack.
}


\examples{
# Create two nested functions
f1 = function() {
    f2();
}
f2 = function() {
    CallStackLevels()
}

f2(); # Returns 1
f1(); # Returns 2

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
