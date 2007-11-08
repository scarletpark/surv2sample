\name{cif}
\alias{cif}
\title{Estimates of Cumulative Incidence Functions for
Competing Risks Data}
\description{
Compute estimates of cumulative incidence functions for one or
more samples of censored data with several competing risks
(types of failure).
}
\usage{
cif(x, group)
}
\arguments{
  \item{x}{an object of class \code{"\link{Survcomp}"}.}
  \item{group}{a vector giving for each observation the group
  to which the observation belongs. If there are \eqn{K} samples,
  \code{group} may only contain values \eqn{1,\dots,K}.
  If \code{group} is missing, observations are assumed to be
  one sample.}
}
\details{
The cumulative incidence function \eqn{F(t,j)} is defined as the probability
of failure by time \eqn{t} from a particular cause \eqn{j} in
the presence of other competing risks \eqn{1,\dots,j-1,j+1,\dots,J},
that is, \eqn{F(t,j)=P(T\leq t, \epsilon=j)}{F(t,j)=P(T<=t, e=j)},
where \eqn{T} is the failure time and \eqn{\epsilon}{e} is the failure
cause.

The function \code{cif} estimates cumulative incidence functions for
all causes and all groups.
}
\value{
\code{cif} returns a list of class \code{"cif"} with the following components:
  \item{[[\var{k}]]}{a list containing estimates for the \eqn{k}th
  sample. It has components \code{time} (sorted times in the
  group \eqn{k}), \code{surv} (the Kaplan--Meier estimate of the
  overall survival in this group), and \code{f}, which is
  a matrix containing estimates of cumulative incidence functions
  for the \eqn{k}th sample. The \eqn{j}th column of \code{f} is
  the cumulative incidence function for the cause \eqn{j}.}
  \item{time}{the vector of sorted times.}
  \item{event}{the vector of corresponding event types.}
  \item{group}{the vector of corresponding group indicators.}
  \item{ncauses}{the number of causes present in the data.}
  \item{ngroups}{the number of samples.}
}
\references{
Klein, J. P. and Moeschberger, M. L. (2003) \emph{Survival Analysis.
Techniques for Censored and Truncated Data.} Springer, New York.
(Section 2.7)
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{\code{\link{Survcomp}}, \code{\link{plot.cif}}}
\examples{
## bone marrow transplant data
data(bmt1)

print(a <- cif(Survcomp(bmt1$time, bmt1$event), bmt1$donor))
str(a)

## several first times and cifs for group 1 (HLA-identical
## sibling donors)
head(cbind(a[[1]]$time, a[[1]]$f))
## several last times and cifs for group 2 (HLA-matched
## unrelated donors)
tail(cbind(a[[2]]$time, a[[2]]$f))

## plot aggregate cumulative incidence functions for each
## donor type to see probabilities within groups
plot(a)
}
\keyword{survival}
