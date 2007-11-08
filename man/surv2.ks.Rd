\name{surv2.ks}
\alias{surv2.ks}
\title{Two-Sample Kolmogorov--Smirnov, Cramer--von Mises and
Anderson--Darling Test for Censored Data}
\description{
Performs the Kolmogorov--Smirnov, Cramer--von Mises and
Anderson--Darling test to compare the distribution of survival
times in two samples of censored data.
}
\usage{
surv2.ks(x, group, process = "w", approx = "lwy", nsim = 2000,
         nsim.plot = 50)
}
\arguments{
  \item{x}{a \code{"Surv"} object, as returned by the \code{\link{Surv}}
  function.}
  \item{group}{a vector indicating to which group each observation belongs.
  May contain values 1 and 2 only.}
  \item{process}{the type of the test process. Possible values are
  \code{"w"} for the difference of Nelson--Aalen estimates
  (asymptotically a Brownian motion, i.e., \strong{W}iener process),
  \code{"b"} for a transformation of this process (asymptotically
  a Brownian \strong{b}ridge).}
  \item{approx}{the method of approximating the distribution of
  test statistics. Possible values are \code{"lwy"} or \code{"mart"}
  for the martingale-based simulation, \code{"perm"} for permutations,
  \code{"boot"} for the bootstrap.}
  \item{nsim}{the number of simulations (martingale simulations,
  permutations or bootstrap samples).}
  \item{nsim.plot}{the number of simulated paths of the test process
  to be returned (for possible plotting). Must be at most \code{nsim}.}
}
\details{
The function implements tests based on functionals of the difference
\eqn{U(t)=\hat A_2(t)-\hat A_1(t)}{U(t)=A_2(t)-A_1(t)}, where
\eqn{\hat A_j(t)}{A_j(t)} are Nelson--Aalen estimates. This process
(properly normalised) is asymptotically a Brownian motion in transformed
time. If \code{process} is \code{"w"}, the supremum (KS) and integral
(CM, AD) test statistics are computed from this process.
If \code{process} is \code{"b"}, the tests are instead based on the
process \eqn{U(t)(1+\hat v(t)/\hat v(\tau))^{-1}}{U(t)(1+v(t)/v(tau))^(-1)},
which is asymptotically a Brownian bridge in transformed time.
}
\value{
A list with class attributes \code{"surv2.int"} and \code{"lwy.test"},
with main components:
  \item{stat.ks}{the Kolmogorov--Smirnov statistic.}
  \item{pval.ks}{the corresponding \eqn{p}-value.}
  \item{pval.ks.asympt}{the asymptotic \eqn{p}-value.}
  \item{stat.cm}{the Cramer--von Mises statistic.}
  \item{pval.cm}{the corresponding \eqn{p}-value.}
  \item{stat.ad}{the Anderson--Darling statistic.}
  \item{pval.ad}{the corresponding \eqn{p}-value.}
  \item{time}{sorted times.}
  \item{test.process}{the test process.}
  \item{test.process.sim}{simulated paths of the test process
  (a matrix with \code{nsim.plot} columns).}
Some of input arguments are also contained in the output.
}
\references{
Andersen, P. K., Borgan, O., Gill, R. D. and Keiding, N. (1993)
\emph{Statistical Models Based on Counting Processes.} Springer, New York. 

Fleming, T. R. and Harrington, D. P. (1991) \emph{Counting Processes
and Survival Analysis.} Wiley, New York.
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{
See the \code{\link[=plot.lwy.test]{plot}} method inherited from
the class \code{"lwy.test"}.

See also \code{\link{surv2.neyman}}, \code{\link{surv2.logrank}},
\code{\link{survdiff}}, \code{\link{survfit}.}
}
\examples{
## gastric cancer data
data(gastric)

## print results
print(a <- surv2.ks(Surv(gastric$time, gastric$event),
    gastric$treatment))
## plot the test process and simulated paths
plot(a)
}
\keyword{survival}
