\name{proprate2.neyman}
\alias{proprate2.neyman}
\alias{summary.proprate2.neyman}
\title{Neyman's Smooth Test of Proportional Rates in Two Samples of Censored Data}
\description{
Checks the assumption of proportional rates (proportional hazards, proportional
odds) in two samples of right-censored data using (possibly data-driven)
Neyman's smooth test.
}
\usage{
proprate2.neyman(x, group, model = 0, data.driven = TRUE,
                 d = ifelse(data.driven, 5, 3),
                 basis = "legendre", time.transf = "F",
                 beta.init = 0, maxiter = 20, eps = 1e-09,
                 choltol = 1e-07)

\method{summary}{proprate2.neyman}(object, ...)
}
\arguments{
  \item{x}{a \code{"Surv"} object, as returned by the \code{\link{Surv}}
  function.}
  \item{group}{a vector indicating to which group each observation belongs.
  May contain values 1 and 2 only.}
  \item{model}{the type of model. Possible values are 0 for proportional
  hazards, 1 for proportional odds.}
  \item{data.driven}{Should the test be data-driven?}
  \item{d}{the number of basis functions for the test with fixed dimension,
  the maximum number of basis functions for the data-driven test.}
  \item{basis}{the basis of functions. Possible values are \code{"legendre"}
  for Legendre polynomials and \code{"cos"} for cosines.}
  \item{time.transf}{the time transformation for basis functions.
  Possible values are \code{"F"} for the 
  distribution function (\eqn{F(t)/F(\tau)}{F(t)/F(tau)}) estimated from the
  pooled sample (recommended), \code{"A"} for the 
  cumulative hazard (\eqn{A(t)/A(\tau)}{A(t)/A(tau)}) and \code{"I"}
  for no transformation (the linear transformation \eqn{t/\tau}{t/tau}).}
  \item{beta.init}{the initial parameter value for iteration in the
  simplified partial likelihood estimation.}
  \item{maxiter}{the maximum number of iterations.}
  \item{eps}{the convergence tolerance parameter. The convergence criterion
  is \eqn{|(l-l_\mathrm{old})/l|<}{|(l-l_old)/l|<}\code{eps}.}
  \item{choltol}{a tolerance parameter for the Cholesky decomposition.}
  \item{object}{an object of class \code{"proprate2.neyman"}, as returned
  by \code{proprate2.neyman}.}
  \item{\dots}{further parameters for printing.}
}
\details{
This function tests the hypothesis that transformation rates
(currently hazard rates or odds functions) are proportional
(their ratio is constant in time) in two samples of censored
survival data.

The proportional rate model is estimated by a two-sample simplification
of the partial likelihood. Then Neyman's smooth test of fit is performed.
In general, Neyman's smooth tests are based on embedding the null hypothesis
in a \eqn{d}-dimensional alternative. Here it consists of expressing the
logarithm of the possibly time-varying ratio of rates as a linear combination
of \eqn{d} basis functions (Legendre polynomials or cosines) in transformed
time. Their significance is tested by a score test. The score is derived
from the simplified partial likelihood. See Kraus (2007a) for details.
The quadratic test statistic is asymptotically chi-square distributed
with \eqn{d} degrees of freedom.

A data-driven choice of the number of basis functions is possible.
The selection is based on a Schwarz-type criterion which is the
maximiser of penalised score statistics for dimensions \eqn{1,\dots,d}.
For the \eqn{p}-value of the data-driven test a two-term approximation
is used (see Kraus (2007b), eq. (12)), as the asymptotic chi-square with
1 d.f. is inaccurate.

If the test is data-driven, the \code{summary} method prints details
on the selection procedure (statistics and penalised statistics for
each dimension). This is equivalent to \code{print(x, detail=TRUE, ...)}.
}
\value{
\code{proprate2.neyman} returns a list of class \code{"proprate2.neyman"}
and \code{"neyman.test"}. Its main components are:
  \item{stat}{the test statistic.}
  \item{pval}{the \eqn{p}-value (based on the chi-square distribution for
  the fixed-dimension test and on the two-term approximation for the
  data-driven test).}
  \item{stats, stats.penal}{statistics and penalised statistics for
  dimensions \eqn{1,\dots,d} (only for data-driven tests).}
  \item{S.dim}{the selected dimension (only for data-driven tests).}
Most input parameters and some further components are included.
}
\references{
Bagdonavicius, V. and Nikulin, M. (2000) On goodness-of-fit
for the linear transformation and frailty models. \emph{Statist.
Probab. Lett.} \bold{47}, 177--188.

Kraus, D. (2007a) Checking proportional rates in the two-sample
transformation model. Research Report 2203, Institute of
Information Theory and Automation, Prague. Available at
\url{http://www.davidkraus.net/surv2sample/}.

Kraus, D. (2007b) Data-driven smooth tests of the proportional hazards
assumption. \emph{Lifetime Data Anal.} \bold{13}, 1--16.
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{
\code{\link{proprate2.ks}}, \code{\link{proprate2.gs}} for other
tests of the proportional rate assumption

\code{\link{proprate2}} for estimation
}
\examples{
## chronic active hepatitis data
data(hepatitis)

## Neyman's test of proportional odds
## test with fixed dimension
proprate2.neyman(Surv(hepatitis$time, hepatitis$status),
    hepatitis$treatment, model = 1, data.driven = FALSE)
## data-driven test
print(a <- proprate2.neyman(Surv(hepatitis$time, hepatitis$status),
    hepatitis$treatment, model = 1, data.driven = TRUE))
## details of the selection procedure
summary(a)
}
\keyword{survival}
