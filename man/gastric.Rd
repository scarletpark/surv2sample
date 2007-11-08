\name{gastric}
\alias{gastric}
\docType{data}
\title{Gastric Cancer Data}
\description{
Survival data from a trial comparing chemotherapy versus combined
chemotherapy plus radiotherapy in the treatment of gastric cancer.
}
\usage{data(gastric)}
\format{
  A data frame with 90 observations (45 in each treatment group)
  with the following 3 variables:
  \describe{
    \item{\code{time}}{survival or censoring time, in days}
    \item{\code{event}}{1 = death, 0 = censoring}
    \item{\code{treatment}}{1 = chemotherapy, 2 = chemotherapy plus radiotherapy}
  }
}
\source{
Stablein, D. M. and Koutrouvelis, I. A. (1985) A two-sample test sensitive
to crossing hazards in uncensored and singly censored data.
\emph{Biometrics} \bold{41}, 643--652. (Page 649)
}
\examples{
data(gastric)
plot(survfit(Surv(time, event) ~ treatment, data = gastric))
}
\keyword{datasets}
