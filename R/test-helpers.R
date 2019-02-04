#' Retrieves an example FacileDataSet
#'
#' A subset of the TCGA data from the BLCA and COAD indications is provided
#' as a FacileDataSet.
#'
#' @export
exampleFacileDataSet <- function() {
  fn <- system.file('extdata', 'exampleFacileDataSet', package='FacileData')
  FacileDataSet(fn)
}

#' Fetches exemplar data for unit testing
#'
#' @export
#' @rdname test-helpers
example_sample_covariates <- function() {
  pdat <- system.file("testdata", "test-sample-covariates.rds",
                      package = "FacileData")
  readRDS(pdat)
}

#' @export
#' @rdname test-helpers
#' @param file.path If `TRUE`, returns the path to the yaml file, otherwise
#'   returns the list-of-list meta definition.
#' @return Either the list-of-list meta definition, or path to the `meta.yaml`
#'   file where these are defined.
example_meta <- function(file.path=FALSE) {
  out <- system.file("testdata", "expected-meta.yaml",
                     package = "FacileData")
  if (!isTRUE(file.path)) {
    out <- yaml::read_yaml(out)
  }
  out
}

#' @export
#' @importFrom yaml read_yaml
#' @rdname test-helpers
#' @return the list-of-list definitions for the example `pData` returned from
#'   [example_sample_covariates()]
example_sample_covariate_definitions <- function() {
  out <- example_meta(file.path=FALSE)
  out$sample_covariate
}
