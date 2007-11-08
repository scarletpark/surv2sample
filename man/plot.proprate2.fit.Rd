\name{plot.proprate2.fit}
\alias{plot.proprate2.fit}
\title{Plot the Two-Sample Proportionl Rate Model for Censored Data}
\description{
This function plots estimates of cumulative rates (cumulative hazards
or odds functions) for two samples of censored data. It may plot both
separate estimates from the two samples and estimates based on the
proportional rate model.
}
\usage{
\method{plot}{proprate2.fit}(x, log.transform = FALSE, diff = FALSE, lwds = 1, cols = 1,
     ltys, ...)
}
\arguments{
  \item{x}{a \code{"proprate2.fit"} object, as returned by \code{\link{proprate2}}.}
  \item{log.transform}{logical. Should the logarithms of cumulative rates be plotted?}
  \item{diff}{logical. Instead of two curves, should the difference of their
  logarithms be plotted?}
  \item{lwds, cols, ltys}{vectors of length equal to the number of curves in plots
  (4 if \code{diff} is \code{FALSE}, 2 if \code{TRUE}). These give line widths,
  colours and line types for each curve. If of length 1, the value is replicated.}
  \item{\dots}{further plotting parameters.}
}
\details{
If \code{diff} is \code{FALSE}, four curves are plotted (two individual sample
estimates and two model based estimates). In this case, \code{ltys} defaults
to \code{c(1,1,2,2)}. If \code{diff} is \code{FALSE}, 
the function plots their differences. Then \code{ltys} defaults to \code{c(1,1)}.
To omit a curve, set the corresponding component of \code{lty} to \code{0}.

Using these plots one may visually assess the validity of the proportional rate
assumption.
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{
\code{\link{proprate2}} for estimation

\code{\link{proprate2.neyman}}, \code{\link{proprate2.ks}},
\code{\link{proprate2.gs}} for tests of the proportional rate
assumption
}
\examples{
## chronic active hepatitis data
data(hepatitis)

## fit the proportional odds model
fit = with(hepatitis, proprate2(Surv(time, status), treatment,
    model = 1))

## plot model-based and model-free estimates of odds functions
plot(fit)
## their logarithms
plot(fit, log.transform = TRUE)
## differences of log-functions
plot(fit, diff = TRUE)
}
\keyword{survival}
