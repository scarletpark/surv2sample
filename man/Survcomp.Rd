\name{Survcomp}
\alias{Survcomp}
\alias{is.Survcomp}
\title{Competing Risks Survival Data Object}
\description{
Create an object representing right-censored data with competing risks
(types of failure).
}
\usage{
Survcomp(time, event)
is.Survcomp(x)
}
\arguments{
  \item{time}{survival times (possibly right-censored).}
  \item{event}{the status indicator, 0 = censored, 1 = dead from cause 1,
  2 = dead from cause 2,\dots Successive integers 1,2,\dots must be used
  for causes.}
  \item{x}{any object.}
}
\details{
This is an extension of the standard function \code{\link{Surv}} of
the package \pkg{survival}. \code{\link{Surv}} allows only death and
censoring, while \code{Survcomp} handles causes of death.
}
\value{
\code{Survcomp} returns a matrix of class \code{"Survcomp"} with two
columns (survival times,
event types) with attribute \code{"ncauses"} containing the number
of different causes of death (which is \code{max(event)}). There is a
print method for objects of class \code{"Survcomp"}, and the function
\code{is.Survcomp} for testing whether an object is of this class.
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{\code{\link{cif}}, \code{\link{Surv}}}
\examples{
## bone marrow transplant data
data(bmt1)

## cause of failure printed in brackets
## censored observations have a "+"
print(a <- Survcomp(bmt1$time, bmt1$event))
is.Survcomp(a)
## what's inside
str(a)
}
\keyword{survival}
