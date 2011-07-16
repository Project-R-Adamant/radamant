\name{3dptpars}
\alias{x.axis3d}
\alias{y.axis3d}
\alias{z.axis3d}
\alias{x.title3d}
\alias{y.title3d}
\alias{z.title3d}
\title{ 3D Plot Axis Formatting}
\alias{getPlotLimits}
\alias{setPlotLimits}

\description{Add and format labels for 3D Plot }
\usage{ 

x.axis3d(xlim = getPlotLimits(1), ylim = getPlotLimits(2), 
zlim = getPlotLimits(3), 
pmat = getProjectionMatrix(), at = NULL,
labels = NULL, theme.params = getCurrentTheme(), 
show.labels = TRUE, grid = theme.params[["xgrid"]],
overrides = list(...), ...) 

y.axis3d(xlim = getPlotLimits(1), ylim = getPlotLimits(2), zlim =
getPlotLimits(3), 
pmat = getProjectionMatrix(), at = NULL,
labels = NULL, theme.params = getCurrentTheme(), 
show.labels = TRUE, grid = theme.params[["ygrid"]],
overrides = list(...), ...) 

z.axis3d(xlim = getPlotLimits(1), ylim = getPlotLimits(2), zlim =
getPlotLimits(3), 
pmat = getProjectionMatrix(), at = NULL, labels =NULL, 
theme.params = getCurrentTheme(), show.labels = TRUE, 
grid = theme.params[["zgrid"]],
overrides = list(...), ...)

x.title3d(xlim = getPlotLimits(1), ylim = getPlotLimits(2), zlim =
getPlotLimits(3), 
pmat = getProjectionMatrix(), title = "", 
theme.params = getCurrentTheme(), ...) 

y.title3d(xlim = getPlotLimits(1), ylim = getPlotLimits(2), zlim =
getPlotLimits(3), 
pmat = getProjectionMatrix(), title = "", 
theme.params = getCurrentTheme(), ...) 

z.title3d(xlim = getPlotLimits(1), ylim = getPlotLimits(2), zlim =
getPlotLimits(3), 
pmat = getProjectionMatrix(), title = "",
theme.params = getCurrentTheme(), ...)

getPlotLimits(which = 1:3, env = getOption("RAdamant"))

setPlotLimits(xlim = NULL
, ylim = NULL
, zlim = NULL
, env = getOption("RAdamant")
)

}
\arguments{
\item{xlim}{ xlim}\item{ylim}{ ylim}\item{zlim}{ zlim}\item{pmat}{
pmat}
\item{at}{at}
\item{which}{which}
\item{env}{environment}
\item{labels}{labels}
\item{title}{ title}
\item{theme.params}{ theme.params}
\item{show.labels}{ show.labels}
\item{grid}{grid}
\item{overrides}{ Overrides list}
\item{...}{Further arguments to or from other methods} 
}
\author{ RAdamant Development Team \email{team@r-adamant.org}}


