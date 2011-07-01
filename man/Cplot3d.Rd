\name{Cplot3d}

\alias{cplot3d}

\title{ 3-Dimensional plotting}

\description{Workhorse function for 3D automatic plotting}

\usage{
cplot3d(x, y, z, fill = c("simple", "colormap", "gradient"), 
main = "", xtitle = "", ytitle = "", ztitle = "", 
xlim = range(x) + 0.1*diff(range(x))*c(-1, 1),
ylim = range(y) + 0.1*diff(range(y))*c(-1, 1), 
zlim = range(z, na.rm = TRUE) + 0.1*diff(range(z, na.rm = TRUE))*c(-1, 1), 
xlabels = NULL, ylabels = NULL,
zlabels = NULL, pre = NULL, post = NULL, 
theme.params = getCurrentTheme(),
overrides = list(...), new.device = FALSE, 
append = FALSE, axis = TRUE,
show.labels = TRUE, show.xlabels = TRUE, 
show.ylabels = TRUE, show.zlabels = TRUE
, ...) 
}

\arguments{ 
\item{x}{ x coordinates for the plot}
\item{y}{ y coordinates for the plot}
\item{z}{ z coordinates for the plot}
\item{fill}{ fill}
\item{main}{main}
\item{xtitle}{ xtitle}
\item{ytitle}{ ytitle}
\item{ztitle}{ ztitle}
\item{xlim}{ xlim}
\item{ylim}{ ylim}
\item{zlim}{ zlim}
\item{xlabels}{ xlables}
\item{ylabels}{ ylabels}
\item{zlabels}{ zlabels}
\item{pre}{ pre}
\item{post}{ post}
\item{theme.params}{ theme.params}
\item{overrides}{ overrides}
\item{new.device}{ new.device}
\item{append}{ append}
\item{axis}{ axis}
\item{show.labels}{ show.labels}
\item{show.xlabels}{ show.xlabels}
\item{show.ylabels}{ show.ylabels}
\item{show.zlabels}{ show.zlabels}
\item{...}{ Further arguments to or from other methods}
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
