\name{cif2.logrank}
\alias{cif2.logrank}
\title{Two-Sample Logrank-Type Test for Cumulative Incidence Functions}
\description{
Compares cumulative incidence functions (CIF) for one failure cause
in two samples of censored competing risks data using Gray's
\eqn{G^\rho}{G^rho}-weighted logrank-type test based on subdistribution
hazards corresponding to CIFs.}
\usage{
cif2.logrank(x, group, cause = 1, rho = 0)
}
\arguments{
  \item{x}{a \code{"Survcomp"} object, as returned by the
  \code{\link{Survcomp}} function.}
  \item{group}{a vector indicating to which group each observation belongs.
  May contain values 1 and 2 only.}
  \item{cause}{For which cause of failure should the CIFs be compared?}
  \item{rho}{the parameter (exponent) of the weight.}
}
\details{
The test compares cumulative incidence functions \eqn{F_1(t,k)},
\eqn{F_2(t,k)} for a particular failure cause \eqn{k}.

This test was proposed by Gray (1988). The statistic is similar
to the weighted logrank test, but instead of ordinary hazards
this test compares subdistribution hazards, which is equivalent
to comparing CIFs (due to the one-to-one relationship between 
subdistribution hazards and CIFs). The \eqn{G^\rho}{G^rho} weight
is of the form \eqn{(1-\hat F_0(t,k))^\rho}{(1-F_0(t,k))^rho},
where \eqn{\hat F_0(t,k)}{F_0(t,k)} is a pooled sample estimate of
the cause \eqn{k}{k} cumulative incidence curve.

This logrank-type test compares subdistribution hazards, and should
not be confused with the ordinary logrank test applied to comparing
cause-specific hazards. Hypotheses on subdistribution hazards (or CIFs)
and cause-specific hazards are not equivalent.
}
\value{
A list of class \code{"cif2.logrank"} with components:
  \item{stat}{the test statistic.}
  \item{pval}{the \eqn{p}-value based on the asymptotic normality.}
Further components are \code{cause} and \code{rho}, the same
as on input.
}
\references{
Gray, R. J. (1988) A class of \eqn{k}-sample tests for comparing the
cumulative incidence of a competing risk. \emph{Ann. Statist.}
\bold{16}, 1141--1154.
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{
\code{\link{cif}} and \code{\link{plot.cif}} for estimation and
plotting of CIFs, \code{\link{cif2.ks}}, \code{\link{cif2.int}} and
\code{\link{cif2.neyman}} for other two-sample tests.
}
\examples{
## bone marrow transplant data
data(bmt1)

## compare CIFs for cause 1 (relapse)
cif2.logrank(Survcomp(bmt1$time, bmt1$event), bmt1$donor, cause = 1)

## compare CIFs for cause 2 (death in remission)
cif2.logrank(Survcomp(bmt1$time, bmt1$event), bmt1$donor, cause = 2)
}
\keyword{survival}
