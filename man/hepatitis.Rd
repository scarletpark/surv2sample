\name{hepatitis}
\alias{hepatitis}
\docType{data}
\title{Chronic Active Hepatitis Data}
\description{
Survival data from a randomised clinical trial of the drug prednisolone
for chronic active hepatitis.
}
\usage{data(hepatitis)}
\format{
A data frame with 44 observations, with the following 3 variables:
  \describe{
    \item{\code{treatment}}{treatment type (1 = drug prednisolone,
    2 = untreated control group)}
    \item{\code{time}}{time following admission to the trial, in months}
    \item{\code{status}}{censoring status (1 = death, 0 = censored)}
  }
}
\source{
Collett, D. (2003) \emph{Modelling Survival Data in Medical Research.
Second Edition.} Chapman \& Hall/CRC. (Appendix D.1, p. 363)\cr
\url{http://www.crcpress.com/e_products/downloads/webdownload/C3251/C3251.zip}
}
\examples{
data(hepatitis)
plot(survfit(Surv(time, status) ~ treatment, data = hepatitis))
}
\keyword{datasets}
