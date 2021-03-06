% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csurv.R
\name{cSurv}
\alias{cSurv}
\alias{as_cSurv}
\alias{as_Surv}
\title{cSurv is a character representation of survival::surv ###}
\usage{
as_cSurv(from)

as_Surv(from)
}
\description{
cSurv serves as a more reliable way to use Surv objects as data.frame columns. A
data.frame is supposed to be able to hold Surv columns. There are multiple special
cases written into base for this. It seems the implementation is incomplete as
subsetting the DF breaks the Surv object. cSurv cannot do anything but get subset
and become a Surv again. In the FacileVerse we hold Surv objects as cSurv, which
allows us to survive a round-trip through an EAV sample metadata table. Survival
analyses can convert cSurv to Surv as needed. It is assumed that all Surv censoring
is right-censored.
}
\examples{
library(survival)
x = Surv(c(14,12,3), event = c(1,0,1))
y = as(x,"cSurv")
z = as(y, "Surv")
x2 = as.character(x)
z2 = as(x2, "Surv")
a = as(x, "cSurv")
b = as(a, "character")
c = as(b, "cSurv")
d = as(c, "Surv")
}
\concept{cSurv}
