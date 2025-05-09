% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize_assay_data.R
\name{normalize_assay_data}
\alias{normalize_assay_data}
\title{Helper functions to normalize assay data into log2 space.}
\usage{
normalize_assay_data(
  x,
  features,
  samples,
  batch = NULL,
  log = TRUE,
  prior.count = 0.1,
  main = NULL,
  verbose = FALSE,
  ...,
  .fds = NULL
)
}
\arguments{
\item{x}{A matrix of raw/unnormalized assay data retrieved from
within the \code{fetch_assay_data()} itself.}

\item{features}{a feature descriptor data.frame that includes the
feature_id's of the rows in \code{x}, as well as the assay name/type they
were pulled from. We assert that all features come from the same assay
type, and the rows here match 1:1 the rows in \code{x}.}

\item{samples}{a sample descriptor for the columns in \code{x}. Rows here
should match columns in \code{x} 1:1.}

\item{batch, main}{parameters sent to \code{\link[=remove_batch_effects]{remove_batch_effects()}} after assay
data has been materialized.}
}
\description{
This is defined for the assay_types defined within this package. If you are
writing a package to handle new types of data, you need to define a
\code{normalize_assay_matrix.ASSAY_TYPE} function. This is experimental.
}
