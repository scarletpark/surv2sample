\name{plot.cif}
\alias{plot.cif}
\title{Plots of Cumulative Incidence Functions for
Competing Risks Data}
\description{
Plot cumulative incidence functions for one or more samples of
censored data with several competing risks (types of failure).
}
\usage{
\method{plot}{cif}(x, by = "group", aggreg.cif = TRUE, orient = "land",
     lwds = 1, cols = 1, ltys = if ((by == "cause") || (by == "c"))
     rep(1:6, len = x$ngroups) else 1, xlab = "", ylab = "",
     ylim, mfrow, mfcol, mains, ...)
}
\arguments{
  \item{x}{an object of class \code{"cif"}, as returned by \code{\link{cif}}.}
  \item{by}{If \code{"group"} or \code{"g"}, for every group one plot containing
  all cumulative incedence functions is produced. If \code{"cause"} or \code{"c"},
  for every cause one plot conatining CIFs for all groups is produced.}
  \item{aggreg.cif}{logical. Used only when curves are plotted by group. If
  \code{aggreg.cif} is \code{TRUE}, a summary plot with aggregate CIFs in
  each group is produced. This means that in each plot the lowest curve
  is the CIF for cause 1, the second lowest is the sum of CIFs for causes 1
  and 2, etc. The area above the top curve is the disease-free survival.}
  \item{orient}{If \code{"land"}, multiple plots are arranged to give
  a landscape plot. Otherwise a portrait orientation is assumed.
  \code{orient} is ignored, if \code{mfrow} or \code{mfcol} is
  present.}
  \item{lwds, cols, ltys}{vectors of length equal to the number of curves in plots (number of groups
  or number of causes). These give line widths, colours and line types for each curve.}
  \item{xlab,ylab}{labels for axes.}
  \item{ylim}{limits on the vertical axis. If missing, determined automatically.}
  \item{mfrow, mfcol}{parameters determining the arrangement of
  plots, passed to \code{par}.}
  \item{mains}{the main title(s) for plots. If missing, titles for
  plots are automatically produced. If \code{mains} is present, it
  must be of length equal to the number of plots or 1 (in which case it
  is replicated).}
  \item{\dots}{other standard parameters for plotting.}
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{\code{\link{cif}}}
\examples{
## bone marrow transplant data
data(bmt1)

## plot aggregate cumulative incidence functions for each
## donor type, i.e., plotting by groups, to see the probability
## structure within groups
plot(cif(Survcomp(bmt1$time, bmt1$event), bmt1$donor))

## plot CIFs by causes, to compare donor type effects on each
## cause of failure
plot(cif(Survcomp(bmt1$time, bmt1$event), bmt1$donor), by = "cause")
}
\keyword{survival}
