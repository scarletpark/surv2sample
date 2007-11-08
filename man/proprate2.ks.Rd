\name{proprate2.ks}
\alias{proprate2.ks}
\title{Kolmogorov--Smirnov Test of Proportional Rates in Two Samples of Censored Data}
\description{
Checks the assumption of proportional rates (proportional hazards, proportional
odds) in two samples of right-censored data using
the Kolmogorov--Smirnov test based on the simplified partial likelihood score
process.
}
\usage{
proprate2.ks(x, group, model = 0, nsim = 2000, nsim.plot = 50,
             beta.init = 0, maxiter = 20, eps = 1e-09)
}
\arguments{
  \item{x}{a \code{"Surv"} object, as returned by the \code{\link{Surv}}
  function.}
  \item{group}{a vector indicating to which group each observation belongs.
  May contain values 1 and 2 only.}
  \item{model}{the type of model. Possible values are 0 for proportional
  hazards, 1 for proportional odds.}
  \item{nsim}{the number of simulations to approximate the \eqn{p}-value.
  Must be positive.}
  \item{nsim.plot}{the number of simulated paths of the test process
  to be returned (for possible plotting). Must be at most \code{nsim}.}
  \item{beta.init}{the initial parameter value for iteration in the
  simplified partial likelihood estimation.}
  \item{maxiter}{the maximum number of iterations.}
  \item{eps}{the convergence tolerance parameter. The convergence criterion
  is \eqn{|(l-l_\mathrm{old})/l|<}{|(l-l_old)/l|<}\code{eps}.}
}
\details{
This function tests the hypothesis that transformation rates
(currently hazard rates or odds functions) are proportional
(their ratio is constant in time) in two samples of censored
survival data.

The proportional rate model is estimated by a two-sample simplification
of the partial likelihood. The test then uses the Kolmogorov--Smirnov supremum
statistic based on the simplified partial likelihood score process.
The \eqn{p}-value is computed using the martingale simulation technique.
}
\value{
A list of class \code{"proprate2.ks"} and \code{"lwy.test"},
with main components:
  \item{stat}{the test statistic.}
  \item{pval.sim}{the simulation based \eqn{p}-value.}
  \item{test.process}{the test process.}
  \item{test.process.sim}{simulated paths of the test process
  (a matrix with \code{nsim.plot} columns).}
  \item{time}{sorted times.}
Some of input parameters and further components are included too.
}
\references{
Bagdonavicius, V. and Nikulin, M. (2000) On goodness-of-fit
for the linear transformation and frailty models. \emph{Statist.
Probab. Lett.} \bold{47}, 177--188.
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{
\code{\link[=plot.lwy.test]{plot}} method inherited from
the class \code{"lwy.test"}

\code{\link{proprate2.neyman}}, \code{\link{proprate2.gs}} for other
tests of the proportional rate assumption

\code{\link{proprate2}} for estimation
}
\examples{
## chronic active hepatitis data
data(hepatitis)

## perform the Komogorov--Smirnov test of proportional odds
a = proprate2.ks(Surv(hepatitis$time, hepatitis$status),
    hepatitis$treatment, model = 1)
a
## plot the test process and simulated paths
plot(a)
}
\keyword{survival}
