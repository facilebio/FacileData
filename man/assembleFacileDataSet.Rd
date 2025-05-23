% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test-helpers.R
\name{assembleFacileDataSet}
\alias{assembleFacileDataSet}
\title{Creata a FacileDataSet from the individual assay lists}
\usage{
assembleFacileDataSet(
  name = "TestMultiModalFacileDataSet",
  path = NULL,
  assays = NULL
)
}
\arguments{
\item{name}{The name of the datasets. Defaults to
\code{"TestMultiModalFacileDataSet"}. Users may create different versions to
test (ie. different assays to include, etc) and may want a more
recognizable name for it.}

\item{path}{Local path to create the FacileDataSet directory. If \code{NULL}
(default) a temporary directory will be created via\code{tempfile()} and will
include \code{name} as the directory's prefix.}

\item{assays}{the nums of the assays to include. If \code{NULL} (default), all
assays will be included.}
}
\description{
Creata a FacileDataSet from the individual assay lists
}
