\name{cif2.ks}
\alias{cif2.ks}
\title{Kolmogorov--Smirnov Two-Sample Test for Cumulative Incidence Functions}
\description{
Compares cumulative incidence functions (CIF) for one failure cause
in two samples of censored competing risks data using the
Kolmogorov--Smirnov-type test.
}
\usage{
cif2.ks(x, group, cause = 1, nsim = 2000, nsim.plot = 50)
}
\arguments{
  \item{x}{a \code{"Survcomp"} object, as returned by the
  \code{\link{Survcomp}} function.}
  \item{group}{a vector indicating to which group each observation belongs.
  May contain values 1 and 2 only.}
  \item{cause}{For which cause of failure should the CIFs be compared?}
  \item{nsim}{the number of simulations to approximate the \eqn{p}-value.
  Must be positive.}
  \item{nsim.plot}{the number of simulated paths of the test process
  to be returned (for possible plotting). Must be at most \code{nsim}.}
}
\details{
The test compares cumulative incidence functions \eqn{F_1(t,k)},
\eqn{F_2(t,k)} for a particular failure cause \eqn{k}.

The test statistic is the maximum absolute difference of the two
cumulative incidence functions. Its asymptotic distribution is
complicated, therefore the martingale-based simulation approximation
is employed. See Lin (1997).
}
\value{
A list with class attributes \code{"cif2.int"} and \code{"lwy.test"},
with components:
  \item{stat}{the test statistic.}
  \item{pval.sim}{the simulation based \eqn{p}-value.}
  \item{test.process}{the test process (difference of the two CIFs).}
  \item{test.process.sim}{simulated paths of the test process
  (a matrix with \code{nsim.plot} columns).}
  \item{time}{sorted times.}
Further components are \code{cause}, \code{nsim}, \code{nsim.plot}, the same
as on input.
}
\references{
Lin, D. Y. (1997) Non-parametric inference for cumulative incidence functions
in competing risks studies. \emph{Stat. Med.} \bold{16}, 901--910.
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{
See the \code{\link[=plot.lwy.test]{plot}} method inherited from
the class \code{"lwy.test"}.

See \code{\link{cif}} and \code{\link{plot.cif}} for estimation and
plotting of CIFs, \code{\link{cif2.int}}, \code{\link{cif2.logrank}} and
\code{\link{cif2.neyman}} for other two-sample tests.
}
\examples{
## bone marrow transplant data
data(bmt1)

## compare CIFs for cause 1 (relapse)
## print results
print(a <- cif2.ks(Survcomp(bmt1$time, bmt1$event), bmt1$donor,
      cause = 1))
## plot the test process and simulated paths
plot(a)

## compare CIFs for cause 2 (death in remission)
## print results
print(a <- cif2.ks(Survcomp(bmt1$time, bmt1$event), bmt1$donor,
    cause = 2))
## plot the test process and simulated paths
plot(a)
}
\keyword{survival}
