\name{proprate2}
\alias{proprate2}
\alias{coef.proprate2.fit}
\title{Fitting the Two-Sample Proportional Rate Transformation Model
for Censored Data}
\description{
\code{proprate2} estimates the two-sample proportional rate transformation model
(proportional hazards, proportional odds) for censored data using the simplified
partial likelihood.
}
\usage{
proprate2(x, group, model = 0, beta.init = 0, maxiter = 20,
          eps = 1e-09)
}
\arguments{
  \item{x}{a \code{"Surv"} object, as returned by the \code{\link{Surv}}
  function.}
  \item{group}{a vector indicating to which group each observation belongs.
  May contain values 1 and 2 only.}
  \item{model}{the type of model. Possible values are 0 for proportional
  hazards, 1 for proportional odds.}
  \item{beta.init}{the initial parameter value for iteration.}
  \item{maxiter}{the maximum number of iterations.}
  \item{eps}{the convergence tolerance parameter. The convergence criterion
  is \eqn{|(l-l_\mathrm{old})/l|<}{|(l-l_old)/l|<}\code{eps}.}
}
\details{
This function fits the proportional rate model for two samples
of censored survival data. Currently two most important models are
implemented: proportional hazards and proportional odds. The
estimation procedure is based on a two-sample simplification of the
partial for the two-sample situation, see Bagdonavicius and Nikulin (2000).
(For proportional hazards, this method is the usual partial likelihood.)
}
\value{
A list of class \code{"proprate2.fit"} with main components:
  \item{beta}{the estimate.}
  \item{var}{its variance.}
  \item{G0}{the cumulative baseline rate (at times \code{time}).}
  \item{time}{sorted times.}
  \item{iter}{the number of iterations used.}
  \item{converged}{logical. Did the iterations (appear to) converge?}
  \item{loglik.init}{the simplified partial likelihood at the initial
  value of the parameter.}
  \item{loglik}{the simplified partial likelihood at the estimate.}
  \item{d11}{the derivative of the score.}
  \item{sigma11}{variance of the score (for proportional hazards \code{sigma11}
  equals \code{d11}).}
  \item{G1, G2}{cumulative transformation rates computed separately
  in the two groups (both of length \eqn{n}, at times \code{time}).}
}
\references{
Bagdonavicius, V. and Nikulin, M. (2000) On goodness-of-fit
for the linear transformation and frailty models. \emph{Statist.
Probab. Lett.} \bold{47}, 177--188.

Kraus, D. (2007) Checking proportional rates in the two-sample
transformation model. Research Report 2203, Institute of
Information Theory and Automation, Prague. Available at
\url{http://www.davidkraus.net/surv2sample/}.
}
\author{David Kraus (\url{http://www.davidkraus.net/})}
\seealso{
There is a \code{\link[=plot.proprate2.fit]{plot}} method for objects
returned by \code{proprate2}.

See \code{\link{proprate2.neyman}}, \code{\link{proprate2.ks}},
\code{\link{proprate2.gs}} for tests of the proportional rate assumption.
}
\examples{
## chronic active hepatitis data
data(hepatitis)

## fit the proportional odds model
fit = with(hepatitis, proprate2(Surv(time, status), treatment,
    model = 1))
fit

## plot model-based and model-free estimates of odds functions
plot(fit)
}
\keyword{survival}
