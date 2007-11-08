\name{cif2.int}
\alias{cif2.int}
\title{Two-Sample Integrated-Difference Test for Cumulative Incidence Functions}
\description{
Compares cumulative incidence functions (CIF) for one failure cause
in two samples of censored competing risks data using the test based
on the integrated difference of CIFs.
}
\usage{
cif2.int(x, group, cause = 1, tau, nsim = 0)
}
\arguments{
  \item{x}{a \code{"Survcomp"} object, as returned by the
  \code{\link{Survcomp}} function.}
  \item{group}{a vector indicating to which group each observation belongs.
  May contain values 1 and 2 only.}
  \item{cause}{For which cause of failure should the CIFs be compared?}
  \item{tau}{the upper limit of the integral in the test statistic.
  If missing, defaults to the maximum of times.}
  \item{nsim}{the number of simulations used to approximate the distribution
  of the test statistic. If 0, no simulations are carried out and the
  asymptotic normal approximation is used.}
}
\details{
The test compares cumulative incidence functions \eqn{F_1(t,k)},
\eqn{F_2(t,k)} for a particular failure cause \eqn{k}.

The method is based on the statistic proposed by Pepe (1991) which
is the integral of \eqn{F_2(t,k)-F_1(t,k)} from 0 to
\eqn{\tau}{tau}. The martingale-based simulation technique and the
variance estimator are described in Bajorunaite and Klein (2007).
}
\value{
A list of class \code{"cif2.int"} with components:
  \item{stat}{the test statistic.}
  \item{pval.asympt}{the \eqn{p}-value based on the asymptotic normality.}
  \item{pval.sim}{the \eqn{p}-value based on simulations (if \code{nsim}>0).}
Further components are \code{cause}, \code{tau}, \code{nsim}, the same
as on input.
}
\references{
Bajorunaite, R. and Klein, J. P. (2007) Two-sample tests of the equality
of two cumulative incidence functions. \emph{Comput. Statist. Data Anal.}
\bold{51}, 4269--4281.

Pepe, M. S. (1991) Inference for events with dependent risks in multiple
endpoint studies. \emph{J. Amer. Statist. Assoc.} \bold{86}, 770--778. 
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{
\code{\link{cif}} and \code{\link{plot.cif}} for estimation and
plotting of CIFs

\code{\link{cif2.ks}}, \code{\link{cif2.logrank}} and
\code{\link{cif2.neyman}} for other two-sample tests
}
\examples{
## bone marrow transplant data
data(bmt1)

## compare CIFs for cause 1 (relapse)
cif2.int(Survcomp(bmt1$time, bmt1$event), bmt1$donor, cause = 1)

## compare CIFs for cause 2 (death in remission)
cif2.int(Survcomp(bmt1$time, bmt1$event), bmt1$donor, cause = 2)
}
\keyword{survival}
