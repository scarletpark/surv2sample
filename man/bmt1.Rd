\name{bmt1}
\alias{bmt1}
\docType{data}
\title{Bone Marrow Transplant Data}
\description{
Data on the treatment of leukaemia by the bone marrow transplantation
from two types of donors (HLA-identical sibling, HLA-matched unrelated),
with two competing causes of failure (relapse, death in remission).
}
\usage{data(bmt1)}
\format{
  A data frame with 1607 observations on the following 3 variables:
  \describe{
    \item{\code{time}}{time from bone marrow transplantation, in months}
    \item{\code{event}}{0 = censoring, 1 = relapse (recurrence of the primary disease),
    2 = death in remission (treatment related mortality)}
    \item{\code{donor}}{1 = HLA-identical sibling donor, 2 = HLA-matched unrelated donor
    (HLA = human leukocyte antigen)}
  }
}
\source{
Klein, J. P. and Andersen, P. K. (2005) Regression modeling of competing
risks data based on pseudovalues of the cumulative incidence function.
\emph{Biometrics} \bold{61}, 223--229.

The original dataset was obtained from the website of 
Biometrics (\url{http://biometrics.tibs.org/}), where it is available
as supplementary material to Klein and Andersen's paper (direct link:
\url{http://biometrics.tibs.org/datasets/031209.txt}).
This is a subset (1607 observations) of the original dataset (1715),
108 cases with mismatched unrelated donors were omitted.
}
\examples{
data(bmt1)

## plot aggregate cumulative incidence functions for each
## donor type
with(bmt1, plot(cif(Survcomp(time, event), donor)))
}
\keyword{datasets}
