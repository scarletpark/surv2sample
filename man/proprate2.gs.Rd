\name{proprate2.gs}
\alias{proprate2.gs}
\title{Gill--Schumacher Test of Proportional Rates in Two Samples of Censored Data}
\description{
Checks the assumption of proportional rates (proportional hazards, proportional
odds) in two samples of right-censored data using the Gill--Schumacher test
based on the comparison of two estimates of the rate ratio.
}
\usage{
proprate2.gs(x, group, model = 0, weight1 = "logrank",
             weight2 = "prentice")
}
\arguments{
  \item{x}{a \code{"Surv"} object, as returned by the \code{\link{Surv}}
  function.}
  \item{group}{a vector indicating to which group each observation belongs.
  May contain values 1 and 2 only.}
  \item{model}{the type of model. Possible values are 0 for proportional
  hazards, 1 for proportional odds.}
  \item{weight1, weight2}{weight functions for the ratio estimates.
  Possible values are \code{"logrank"}, \code{"prentice"},
  \code{"gehan"}.}
}
\details{
This function performs the Gill--Schumacher test of the hypothesis
that transformation rates (currently hazard rates or odds functions)
are proportional (their ratio is constant in time) in two samples of
censored survival data.

The test was proposed by Gill and Schumacher (1987) for proportional
hazards, see Kraus (2007) for its extension to proportional transformation
rates. The test statistic compares two weighted estimates of the ratio of
rates. Possible weights are of the logrank, Prentice--Wilcoxon or
Gehan type.
}
\value{
A \code{"proprate2.gs"} object with components:
  \item{stat}{the test statistic.}
  \item{pval}{the \eqn{p}-value.}
  \item{eta1, eta2}{weighted estimates of the rate ratio.}
Some of input parameters are included too.
}
\references{
Gill, R. and Schumacher, M. (1987) A simple test of the proportional
hazards assumption. \emph{Biometrika} \bold{74}, 289--300.

Kraus, D. (2007) Checking proportional rates in the two-sample
transformation model. Research Report 2203, Institute of
Information Theory and Automation, Prague. Available at
\url{http://www.davidkraus.net/surv2sample/}.
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{
\code{\link{proprate2.neyman}}, \code{\link{proprate2.ks}} for other
tests of the proportional rate assumption

\code{\link{proprate2}} for the simplified partial likelihood estimation
}
\examples{
## chronic active hepatitis data
data(hepatitis)

## perform the Gill--Schumacher test of proportional odds
proprate2.gs(Surv(hepatitis$time, hepatitis$status),
    hepatitis$treatment, model = 1)
}
\keyword{survival}
