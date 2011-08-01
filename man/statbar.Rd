\name{statbar}
\alias{statusbar}
\title{Status bar}

\description{Interactive status bar for updating completion percentage to console.}

\usage{
statusbar(message = "Computing.."
	, status = 0
	, n = 1
	, N = 1
	, step = 0.01) 
	}

\arguments{
\item{message}{The message to be sent to console.}
\item{status}{The percentage of completion (status in [0, 1]).}
\item{n}{The current value of the loop counter.}
\item{N}{The total number of iterations}
\item{step}{The percentage increment by which the status is updated.}
}

\value{The updated status for the next iteration.}

\details{This function is meant to be used inside a loop, to inform the user about the current status of the processing.}

\examples{
# Number of iterations
N = 1000;
# Set the message
msg = "Still running..";
# Init Status bar
status = 0;
# Set the step to 0.05. The status bar is updated by 5% each time
step = 0.05;

# Start looping
for(n in 1:N) {
    # Do something
    # ... some code ...
    
    # Update the status (note how status is reused at each iteration)
    status = statusbar(message = msg
                            , status = status
                            , n = n
                            , N = N
                            , step = step);
}

}
\author{RAdamant Development Team \email{team@r-adamant.org}}



