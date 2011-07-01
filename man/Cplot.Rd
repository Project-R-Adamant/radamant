\name{Cplot}

\alias{cplot}

\title{ 2-Dimensional Plotting}
\description{Workhorse function for automatic plotting}

\usage{cplot(X, base = NULL, xrange = NULL, yrange = NULL, theme.params =
getCurrentTheme(), xtitle = "" , xlabels = NULL , ytitle = "" , ylabels =
NULL	, ytitle2 = "" , ylabels2 = NULL, show.xlabels = TRUE , show.ylabels =
TRUE , main = "" , legend = NULL , legend.col = theme.params[["col"]],
show.legend = TRUE , shaded = FALSE, grid = TRUE , overrides = list(...),
new.device = FALSE,append = FALSE, multicolor = FALSE , ...) }

\arguments{ 
\item{X}{ x coordinates of the plot}
\item{base}{ y coordinates of the plot}
\item{xrange}{ x axis range}
\item{yrange}{ y axis range}
\item{theme.params}{ Retrieve RAdamant graphical theme}
\item{xtitle}{ title for the x axis}
\item{xlabels}{ labels for x tick marks}
\item{ytitle}{ title for the y axis}
\item{ylabels}{ labels for y tick marks}
\item{ytitle2}{ title for the right-y axis}
\item{ylabels2}{ labels for right-y tick marks}
\item{show.xlabels}{ LOGICAL. Show labels on the x axis}
\item{show.ylabels}{ LOGICAL. Show labels on the y axis}
\item{main}{ Main title for the plot}
\item{legend}{ Add the legend}
\item{legend.col}{ Colors for the elements in the legend}
\item{show.legend}{ LOGICAL. Display the legend in the plot}
\item{shaded}{ LOGICAL. Insert shaded under the plot}
\item{grid}{ LOGICAL. Draw a grid}
\item{overrides}{ overrides list}
\item{new.device}{ open new.device window}
\item{append}{ append to existing plot}
\item{multicolor}{ multiple colors}
\item{...}{ additional arguments for generic funciotn "plot"}
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
