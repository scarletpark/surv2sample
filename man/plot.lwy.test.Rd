\name{plot.lwy.test}
\alias{plot.lwy.test}
\title{Plot an Observed Test Process and Simulated Paths}
\description{
This function plots an observed test process along with simulated
paths. Such a plot may be used to visually assess the hypothesis.
}
\usage{
\method{plot}{lwy.test}(x, lwds = c(3, 1), cols = c("black", "grey"),
     ltys = c(1, 1), ...)
}
\arguments{
  \item{x}{an object of class \code{"lwy.test"}.}
  \item{lwds}{a vector of length 2, giving the line width for the
  observed process and for simulated processes.}
  \item{cols}{a vector of length 2, giving the line colour for the
  observed process and for simulated processes.}
  \item{ltys}{a vector of length 2, giving the line type for the
  observed process and for simulated processes.}
  \item{\dots}{other standard parameters for plotting.}
}
\details{
In survival analysis, \eqn{p}-values of Kolmogorov--Smirnov and other tests
based on a test process may be approximated by the martingale-based
simulation technique originally proposed by Lin, Wei and Ying (1993).
By plotting the observed test process along with a suitable number of
simulated paths, one may visually assess the validity of the hypothesis.

This function is a general plot method for such tests. Results
of such tests are objects with class attribute \code{"lwy.test"},
at least containing sorted times in \code{x$time}, the observed process
in \code{x$test.process} and a number (\code{x$nsim.plot}) of simulated
processes in \code{x$test.process.sim}. Results of functions listed in
\sQuote{See Also} below are handled by this plot method.
}
\references{
Lin, D. Y., Wei, L. J. and Ying, Z. (1993) Checking the Cox model with
cumulative sums of martingale-based residuals. \emph{Biometrika} \bold{80},
557--572.
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{\code{\link{surv2.ks}}, \code{\link{cif2.ks}}, \code{\link{proprate2.ks}}}
\examples{
## gastric cancer data
data(gastric)

## Kolmogorov--Smirnov test of equal survival distributions
## test process = difference of Nelson--Aalen estimates
## plot the observed test process and simulated paths
plot(surv2.ks(Surv(gastric$time, gastric$event),
    gastric$treatment, approx="mart"))
}
\keyword{survival}
