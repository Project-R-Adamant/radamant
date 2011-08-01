\name{splitwdw}
\alias{splitWindow}
\title{Split Window}
\description{Given an input size N, splits the sequence 1:N into sliding or extended windows and returnn the endpoint indexes of each window.}

\usage{
splitWindow(N
	, direction = c("forward", "backward")
	, mode = c("EW", "SW")
	, from = NULL
	, win.size = 1
	, by = 1
	, labels = 1:N
	, ...
	) 
}
\arguments{
\item{N}{The size of the entire window to be split}
\item{direction}{Controls on which direction the next sub-window is computed. One of "forward" or "backward".}
\item{mode}{Controls how windows endpoint indexes are computed. 
If "EW" (Extended Windows), starting with an initial window of size win.size 
at each step the previous sub-window is extended with additional 'by' points on the side specified by 'direction'.\cr
If "SW" (Sliding Windows), the size on the windows is constant: at each step the previous sub-window is shifted on by the quantity 'by' on the side specified by 'direction'.
}
\item{from}{The starting point from wich the first window is calculated}
\item{win.size}{The initial size of the first window if mode = "EW". The size of all windows if mode = "SW"}
\item{by}{Controls the amount of extension or shift (depending on the mode parameters) of the windows.}
\item{labels}{The labels associated to the N data points of the full window.}
\item{...}{Further arguments to or from other methods.}
}

\value{
A matrix with columns [start.idx, end.idx]. Each row represents the endpoints indexes of a corresponding sub-window.
}

\examples{
## Forward Extended Window
splitWindow(N = 30
	# Start with a window of size 3
	, win.size = 3
	# Start from position 5
	, from = 10
	# Move forward
	, direction = "forward"
	# Extended mode
	, mode = "EW"
	# Increase the size by 5 at each step
	, by = 5
	)

## Backward Extended Window
splitWindow(N = 30
	# Start with a window of size 3
	, win.size = 3
	# Start from position 20
	, from = 20
	# Move backward
	, direction = "backward"
	# Extended mode
	, mode = "EW"
	# Increase the size by 2 at each step
	, by = 2
	)

## Forward Sliding Window
splitWindow(N = 30
	# windows of size 5
	, win.size = 5
	# Move forward
	, direction = "forward"
	# Sliding mode
	, mode = "SW"
	# Slide forward by 5 at each step. This produces non overlapping windows.
	, by = 5
	)

## Backward Sliding Window
splitWindow(N = 30
	# windows of size 3
	, win.size = 3
	# Move backward
	, direction = "backward"
	# Sliding mode
	, mode = "SW"
	# Slide backward by 5 at each step. 
	, by = 5
	)

}

\author{RAdamant Development Team \email{team@r-adamant.org}}



