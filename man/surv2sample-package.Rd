\name{surv2sample-package}
\alias{surv2sample-package}
\alias{surv2sample}
\docType{package}
\title{
Two-Sample Tests for Survival Analysis
}
\description{
The package provides tests for comparing two survival
  distributions, testing equality of two cumulative incidence
  functions under competing risks and checking goodness of fit of
  proportional rate models (proportional hazards, proportional
  odds) for two samples.

}
\details{
%\tabular{ll}{
%Package: \tab surv2sample\cr
%Type: \tab Package\cr
%Version: \tab 0.1-1\cr
%Date: \tab 2007-11-15\cr
%License: \tab GPL (>= 2)\cr
%}
The package implements many two-sample tests for right-censored
survival data. Three main areas and corresponding methods are:
\itemize{
\item comparison of two survival distributions
\itemize{
\item \code{\link{surv2.logrank}}: weighted logrank tests and their combinations (max, sum)
\item \code{\link{surv2.neyman}}: Neyman's smooth test and its data-driven version
\item \code{\link{surv2.ks}}: Kolmogorov--Smirnov, Cramer--von Mises and
Anderson--Darling test
}
\item comparison of two cumulative incidence functions for competing
risks data
\itemize{
\item \code{\link{cif}}: estimation and plotting of cumulative incidence functions
\item \code{\link{cif2.logrank}}: logrank-type test for subdistribution hazards
\item \code{\link{cif2.neyman}}: Neyman's smooth test and its data-driven version
\item \code{\link{cif2.ks}}: Kolmogorov--Smirnov test
\item \code{\link{cif2.int}}: integrated-difference test
}
\item goodness of fit tests of the proportional rate assumption
(proportional hazards or proportional odds functions in two samples)
\itemize{
\item \code{\link{proprate2}}: estimation based on the simplified partial likelihood
\item \code{\link{proprate2.ks}}: Kolmogorov--Smirnov test
\item \code{\link{proprate2.neyman}}: Neyman's smooth test and its data-driven version
\item \code{\link{proprate2.gs}}: Gill--Schumacher type test
}
}

The package also provides three datasets:
\itemize{
\item \code{\link{gastric}}: gastric cancer data
\item \code{\link{bmt1}}: bone marrow transplant data
\item \code{\link{hepatitis}}: chronic active hepatitis data
}
}
\author{
David Kraus (\url{http://www.davidkraus.net/})
}
\references{
There are several technical reports on the topics covered by this package
available at \url{http://www.davidkraus.net/surv2sample/}.
}
\keyword{package}
\keyword{survival}
