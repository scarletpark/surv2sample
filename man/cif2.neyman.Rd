\name{cif2.neyman}
\alias{cif2.neyman}
\alias{summary.cif2.neyman}
\title{Two-Sample Neyman's Smooth Test for Cumulative Incidence Functions}
\description{
Compares cumulative incidence functions (CIF) for one failure cause
in two samples of censored competing risks data using (possibly
data-driven) Neyman's smooth test.
}
\usage{
cif2.neyman(x, group, cause = 1, data.driven = FALSE,
            d = ifelse(data.driven, 5, 3), basis = "legendre",
            choltol = 1e-07)

\method{summary}{cif2.neyman}(object, ...)
}
\arguments{
  \item{x}{a \code{"Survcomp"} object, as returned by the
  \code{\link{Survcomp}} function.}
  \item{group}{a vector indicating to which group each observation belongs.
  May contain values 1 and 2 only.}
  \item{cause}{For which cause of failure should the CIFs be compared?}
  \item{data.driven}{Should the test be data-driven?}
  \item{d}{the number of basis functions for the test with fixed dimension,
  the maximum number of basis functions for the data-driven test.}
  \item{basis}{the basis of functions. Possible values are \code{"legendre"}
  for Legendre polynomials and \code{"cos"} for cosines.}
  \item{choltol}{a tolerance parameter for the Cholesky decomposition.}
  \item{object}{an object of class \code{"cif2.neyman"}, as returned
  by the function \code{cif2.neyman}.}
  \item{\dots}{further parameters for printing.}
}
\details{
The test compares cumulative incidence functions \eqn{F_1(t,k)},
\eqn{F_2(t,k)} for a particular failure cause \eqn{k}.

Neyman-type smooth tests are based on embedding the null hypothesis
in a \eqn{d}-dimensional alternative. The embedding is here formulated
in terms of subdistribution hazards derived from CIFs. The log-ratio
of subdistribution hazards is expressed as a combination of \eqn{d}
basis functions (Legendre polynomials or cosines) in transformed time,
and their significance is tested by a score test. See Kraus (2007b) for
details. The quadratic test statistic is asymptotically chi-square
distributed with \eqn{d} degrees of freedom.

A data-driven choice of the number of basis functions is possible.
The selection is based on a Schwarz-type criterion which is the
maximiser of penalised score statistics for dimensions \eqn{1,\dots,d}.
For the \eqn{p}-value of the data-driven test a two-term approximation
is used (see Kraus (2007a), eq. (12)), as the asymptotic chi-square with
1 d.f. is inaccurate.

If the test is data-driven, the \code{summary} method prints details
on the selection procedure (statistics and penalised statistics for
each dimension). This is equivalent to \code{print(x, detail=TRUE, ...)}.
}
\value{
\code{cif2.neyman} returns a list with class attributes \code{"cif2.neyman"}
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
Kraus, D. (2007a) Data-driven smooth tests of the proportional hazards
assumption. \emph{Lifetime Data Anal.} \bold{13}, 1--16.

Kraus, D. (2007b) Smooth tests of equality of cumulative incidence functions
in two samples. Research Report 2197, Institute of Information Theory and
Automation, Prague. Available at \url{http://www.davidkraus.net/surv2sample/}.
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{
\code{\link{cif}} and \code{\link{plot.cif}} for estimation and
plotting of CIFs, \code{\link{cif2.ks}}, \code{\link{cif2.int}} and
\code{\link{cif2.logrank}} for other two-sample tests.
}
\examples{
## bone marrow transplant data
data(bmt1)

## compare CIFs for cause 1 (relapse)
## test with fixed dimension
cif2.neyman(Survcomp(bmt1$time, bmt1$event), bmt1$donor, cause = 1,
    data.driven = FALSE)
## data-driven test
print(a <- cif2.neyman(Survcomp(bmt1$time, bmt1$event), bmt1$donor,
    cause = 1, data.driven = TRUE))
## print details on the selection procedure
summary(a)

## compare CIFs for cause 2 (death in remission)
cif2.neyman(Survcomp(bmt1$time, bmt1$event), bmt1$donor, cause = 2)
}
\keyword{survival}
