% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api.R
\name{has_assay}
\alias{has_assay}
\title{Checks facile objects to see if they have a particular assay name.}
\usage{
has_assay(x, assay_name = NULL, ...)
}
\description{
If \code{x} a FacileDataStore, return TRUE/FALSE if the assay is stored inside
If \code{x} is a facile_frame, it will append a \verb{has_<assay_name>} with
\code{TRUE}/\code{FALSE} indicating wether or not a particular sample has assay data.
}
