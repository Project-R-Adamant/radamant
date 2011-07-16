\name{mdbuglev}
\alias{getDebugLevel}
\alias{setDebugLevel}
\title{Manage Debug Level}
\description{
Set and retrieve the level of debugging. 
Control how much information is sent to the log about the execution of each function executed.
}

\usage{
setDebugLevel(level = 1, env = getOption("RAdamant"))
getDebugLevel(env = getOption("RAdamant"))
}

\arguments{
\item{level}{The level of debug required (level >= 0). See details.}
\item{env}{The environment where the info is stored (DEFAULT = getOption("RAdamant")).}
}

\value{
The current level of debugging.
}

\details{
The amount of information sent to log depends on the debug level:
\itemize{
\item{level = 0: No information is sent to the log.}
\item{level = 1: Information about main body and conditional executions.}
\item{level = 2: Include information about first level inner loop.}
\item{level = 3: Include information about second level inner loop (loop within loop).}
\item{level = N: Include information about N-th level inner loop.}
}
}

\examples{
# Retrieve current debug level
getDebugLevel();

# Set minimal level of debugging and traceback
setDebugLevel(1);
setDebugTraceLevel(1);
# Enable Console logging
setConsoleLogging(TRUE);

# Compute FFT on some random two-colums matrix. Prints nothing because FFT.default has no logging message
x = FFT(matrix(cumsum(rnorm(256)), 128, 2), plot = FALSE)
plot(x, shaded = FALSE) # Prints nothing because plot.default has no logging message

# Increase Traceback level
setDebugTraceLevel(2);
# Now prints logging info for plot.FFT
plot(x, shaded = FALSE) 

# Increase Debug level
setDebugLevel(2);
# Now prints additional logging info for plot.FFT (from code executed inside a loop)
plot(x, shaded = FALSE) 

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
