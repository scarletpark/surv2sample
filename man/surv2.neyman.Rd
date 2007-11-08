\name{surv2.neyman}
\alias{surv2.neyman}
\alias{summary.surv2.neyman}
\title{Two-Sample Neyman's Smooth Test for Censored Data}
\description{
Compares survival distributions in two samples of censored data
using (possibly data-driven) Neyman's smooth test.
}
\usage{
surv2.neyman(x, group, data.driven = FALSE, subsets = "nested",
             d = ifelse(data.driven, 5, 3), d0 = 0,
             basis = "legendre", time.transf = "F",
             approx = "perm", nsim = 2000, choltol = 1e-07)

\method{summary}{surv2.neyman}(object, ...)
}
\arguments{
  \item{x}{a \code{"Surv"} object, as returned by the \code{\link{Surv}}
  function.}
  \item{group}{a vector indicating to which group each observation belongs.
  May contain values 1 and 2 only.}
  \item{data.driven}{Should the test be data-driven?}
  \item{subsets}{the class of subsets of basis functions amonng which
  the data-driven test selects. Possible values are \code{"nested"}
  and \code{"all"}.}
  \item{d}{the number of basis functions for the test with fixed dimension,
  the maximum number of basis functions for the data-driven test.}
  \item{d0}{the number of high-priority functions for the data-driven test.
  The selection rule selects among subsets containing basis functions
  1,\dots,\code{d0}. For nested subsets, \code{d0} equal to 0 or 1
  is equvalent. For all subsets, \code{d0} equal to 0 means that there
  is no high-priority function and any nonempty subset may be selected.}
  \item{basis}{the basis of functions. Possible values are \code{"legendre"}
  for Legendre polynomials and \code{"cos"} for cosines.}
  \item{time.transf}{the time transformation for basis functions.
  Possible values are \code{"F"} for the 
  distribution function (\eqn{F(t)/F(\tau)}{F(t)/F(tau)})
  (recommended), \code{"A"} for the 
  cumulative hazard (\eqn{A(t)/A(\tau)}{A(t)/A(tau)}) and \code{"I"}
  for no transformation (the linear transformation \eqn{t/\tau}{t/tau}).}
  \item{approx}{the method of approximating the distribution of the
  test statistic. Possible values are \code{"perm"} for permutations,
  \code{"boot"} for the bootstrap, \code{"asympt"} for asymptotics.}
  \item{nsim}{the number of simulations. This means the number of
  permutations or bootstrap samples when \code{approx} is
  \code{"perm"} or \code{"boot"}. When \code{approx} is \code{"asympt"},
  \code{nsim} is the number of simulations to approximate the asymptotic
  distribution (only needed for the data-driven test with all subsets
  and \code{d0} equal to 0).}
  \item{choltol}{a tolerance parameter for the Cholesky decomposition.}
  \item{object}{an object of class \code{"surv2.neyman"}, as returned
  by the function \code{surv2.neyman}.}
  \item{\dots}{further parameters for printing.}
}
\details{
In general, Neyman's smooth tests are based on embedding the null hypothesis
in a \eqn{d}-dimensional alternative. The embedding is here formulated
in terms of hazard functions. The logarithm of the hazard ratio is
expressed as a combination of \eqn{d} basis functions (Legendre polynomials
or cosines) in transformed time, and their significance is tested by a score
test. See Kraus (2007a) for details.
The quadratic test statistic is asymptotically chi-square
distributed with \eqn{d} degrees of freedom.

A data-driven choice of basis functions is possible. The selection is
based on a Schwarz-type criterion which is the maximiser of penalised score
statistics for over a class of nonempty subsets of
\eqn{\{1,\dots,d\}}{{1,\dots,d}}. Either nested subsets with increasing
dimension or all subsets may be used. By choosing \eqn{d_0>0}{d0>0} one
requires that functions with indexes \eqn{\{1,\dots,d_0\}}{{1,\dots,d0}}
be always included in subsets.

If all subsets are used with \eqn{d_0=0}{d0=0}, the data-driven test
statistic is asymptotically distributed as the maximum of (generally
dependent) chi-square variables with 1 d.f. This asymptotic approximation
is accurate. In other cases, the statistic is asymptotically chi-square
distributed with \eqn{d^*=\max(1,d_0)}{d^*=max(1,d0)} degrees of freedom.
For nested subsets with \eqn{d^*=1} a two-term approximation may be used
(see Kraus (2007b), eq. (12)). Otherwise the asymptotics is unreliable.

In any case, one may use permutations or the bootstrap.

If the test is data-driven, the \code{summary} method prints details
on the selection procedure (statistics and penalised statistics for
each considered subset). This is equivalent to \code{print(x, detail=TRUE, ...)}.
}
\value{
A list of class \code{"surv2.neyman"} and \code{"neyman.test"},
with main components:
  \item{stat}{the test statistic.}
  \item{pval}{the \eqn{p}-value.}
  \item{stats, stats.penal}{the score statistic and penalised score
  statistic for each considered subset (only for data-driven tests).}
  \item{S.dim}{the dimension of the selected set (only for data-driven tests).}
  \item{S.set}{the selected set (only for data-driven tests).}
Most input parameters and some further components are included.
}
\references{
Kraus, D. (2007a) Adaptive Neyman's smooth tests of homogeneity
of two samples of survival data. Research Report 2187,
Institute of Information Theory and Automation, Prague.
Available at \url{http://www.davidkraus.net/surv2sample/}.

Kraus, D. (2007b) Data-driven smooth tests of the proportional hazards
assumption. \emph{Lifetime Data Anal.} \bold{13}, 1--16.
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{
\code{\link{surv2.logrank}}, \code{\link{surv2.ks}},
\code{\link{survdiff}}, \code{\link{survfit}}
}
\examples{
## gastric cancer data
data(gastric)

## test with fixed dimension
surv2.neyman(Surv(gastric$time, gastric$event), gastric$treatment,
    data.driven = FALSE)

## data-driven test with nested subsets
## without minimum dimension (i.e., minimum dimension 1)
summary(surv2.neyman(Surv(gastric$time, gastric$event),
    gastric$treatment, data.driven = TRUE, subsets = "nested"))
## with minimum dimension 3
summary(surv2.neyman(Surv(gastric$time, gastric$event),
    gastric$treatment, data.driven = TRUE, subsets = "nested",
    d0 = 3))

## data-driven test with all subsets
## without high-priority functions
summary(surv2.neyman(Surv(gastric$time, gastric$event),
    gastric$treatment, data.driven = TRUE, subsets = "all"))
## with 2 high-priority functions
summary(surv2.neyman(Surv(gastric$time, gastric$event),
    gastric$treatment, data.driven = TRUE, subsets = "all",
    d0 = 2))
}
\keyword{survival}
