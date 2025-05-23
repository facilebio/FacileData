% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample-covariate-summary.R
\name{sample_covariate_summary}
\alias{sample_covariate_summary}
\title{Provides a tibble summary of the covariates available over a sample space.}
\usage{
sample_covariate_summary(x, ..., detailed = FALSE)
}
\arguments{
\item{x}{an object (samples) to summarize covariates over}

\item{detailed}{When \code{TRUE}, summarizes the individual levels of covariates.
Default is \code{FALSE}.}
}
\description{
This is similar to the \verb{FacileData::summary.eav_covariates function}, but
that function is downstream of a \code{collect()}. This function takes in a samples
facile_frame (\strong{preferably a lazy_tbl/tbl_sql}) and summaries the covariates
availalbe there. Keeping this function in SQL space as long as possible
assists in the initial load and presentation of sample filters of monstrously
sized datasets.
}
