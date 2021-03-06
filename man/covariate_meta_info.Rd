% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entity-attribute-value.R
\name{covariate_meta_info}
\alias{covariate_meta_info}
\title{Retrieve the meta information about a covariate for EAV decoding}
\usage{
covariate_meta_info(covariate, .fds, covdefs = NULL)
}
\arguments{
\item{covariate}{the name of the covariate}

\item{.fds}{the \code{FacileDataSet}}

\item{covdefs}{The \code{covariate_definitions(.fds)} list}
}
\value{
a list of covariate information with the following elements:
\verb{$name}, \verb{$type}, \verb{$class}, \verb{$description},
\verb{$label}, \verb{$is.factor}, (and maybe \verb{$levels})
}
\description{
Mappings that define attribute-value encodings into R-native objects are
stored in a \code{FacileDataSet}'s \code{meta.yaml} file, in the \code{sample_covariate}
section.
}
