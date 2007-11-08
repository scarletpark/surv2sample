\name{surv2.logrank}
\alias{surv2.logrank}
\title{Two-Sample Weighted Logrank Tests and Their Combinations}
\description{
Compares the distribution of survival times in two samples of censored
data using the weighted logrank test with the
\eqn{G^{\rho,\gamma}}{G(rho,gamma)} class of weights, or combinations
(maximum or weighted sum) of several weighted logrank statistics.
}
\usage{
surv2.logrank(x, group, rho.gamma = c(0, 0), comb, sum.weights,
              approx = "perm", nsim = 2000, choltol = 1e-07)
}
\arguments{
  \item{x}{a \code{"Surv"} object, as returned by the \code{\link{Surv}}
  function.}
  \item{group}{a vector indicating to which group each observation belongs.
  May contain values 1 and 2 only.}
  \item{rho.gamma}{parameters \eqn{(\rho,\gamma)}{(rho,gamma)} of the weight.
  This must be a vector of length 2 or a list of vectors of length 2.}
  \item{comb}{the method of combining test statistics. Possible values
  are \code{"max"} and \code{"sum"}. \code{comb} is ignored, if there is
  only one set of parameters in \code{rho.gamma}, and defaults to \code{"max"}
  otherwise.}
  \item{sum.weights}{weights for the weighted sum in the \code{"sum"}
  combination method. Defaults to equal weights for all statistics.}
  \item{approx}{the method of approximating the distribution of the
  test statistic. Possible values are \code{"perm"} for permutations,
  \code{"boot"} for the bootstrap, \code{"asympt"} for asymptotics.}
  \item{nsim}{the number of simulations. This means the number of
  permutations or bootstrap samples when \code{approx} is
  \code{"perm"} or \code{"boot"}. When \code{approx} is \code{"asympt"},
  \code{nsim} is the number of simulations to approximate the asymptotic
  distribution (only needed for combination methods).}
  \item{choltol}{a tolerance parameter for the Cholesky decomposition.}
}
\details{
The logrank test uses \eqn{G^{\rho,\gamma}}{G(rho,gamma)} weights defined as
\eqn{\hat S(t)^\rho(1-\hat S(t))^\gamma}{S(t)^rho (1-S(t))^gamma}, where
\eqn{\hat S(t)}{S(t)} is the Kaplan--Meier estimate computed from the pooled sample.

Combination tests are based on a cluster of several weighted logrank
statistics with different parameters \eqn{(\rho,\gamma)}{(rho,gamma)}.
The maximum statistic uses the
maximum of absolute values of standardised statistics. The sum
statistic uses the weighted sum of absolute values of standardised
statistics.
See Chapter 7 of Fleming and Harrington (1991) for details.
}
\value{
A list of class \code{"surv2.logrank"} with components:
  \item{stat}{the test statistic.}
  \item{pval}{the \eqn{p}-value.}
  \item{rho.gamma}{as on input (for a test using a single statistic).}
  \item{stats}{a matrix containing parameters of weights, corresponding
  statistics and \eqn{p}-values (for a combination test).}
Further components include \code{approx}, \code{nsim}, \code{comb},
\code{sum.weights} from input.
}
\references{
Andersen, P. K., Borgan, O., Gill, R. D. and Keiding, N. (1993)
\emph{Statistical Models Based on Counting Processes.} Springer, New York. 

Fleming, T. R. and Harrington, D. P. (1991) \emph{Counting Processes
and Survival Analysis.} Wiley, New York.
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{
\code{\link{surv2.neyman}}, \code{\link{surv2.ks}},
\code{\link{survdiff}}, \code{\link{survfit}}
}
\examples{
## gastric cancer data
data(gastric)

## Prentice--Wilcoxon (G^1) test
surv2.logrank(Surv(gastric$time, gastric$event),
    gastric$treatment, rho.gamma = c(1,0))

## combination of G(0,0), G(1,0), G(1,1) statistics
## maximum test
print(a <- surv2.logrank(Surv(gastric$time, gastric$event),
    gastric$treatment, rho.gamma = list(c(0,0), c(1,0), c(1,1)),
    comb = "max"))
## print individual statistics
a$stats
}
\keyword{survival}
