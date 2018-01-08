#' Creates a FacileDb connection to a test database
#'
#' TODO: This function currently uses the Atezo database as a test. We need to
#'       create a REAL testing database.
#'
#' @export
#' @param db.fn The path on the filesystem to the facile database
#' @param covdef The path on the filesystem to the yaml file that provides
#'   the covariate definitions
exampleFacileDataSet <- function(db.type=c('sqlite', 'monetdblite')) {
  fn <- system.file('extdata', 'exampleFacileDataSet', package='FacileData')
  FacileDataSet(fn)
}

## This uses GenomicsTools to create a test database
createExampleFacileDataSet <- function(out.fn) {
  ## I am farming this out to the FacileTCGA package and creating a smaller
  ## version of the TCGA -- this is cicrular, I know.
  ##
  ## See FacileTCGA/inst/build-test-db
  stop("This is created in FacileTCGA/inst/build-test-db")
}

#' Fetches exemplar pData
#' @export
#' @rdname test-helpers
example_sample_covariates <- function() {
  pdat <- system.file("testdata", "test-sample-covariates.rds",
                      package = "FacileData")
  readRDS(pdat)
}

#' Fetches the meta file for the example FacileDataSet
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

#' Fetches the EAV definitions for the sample covariates
#'
#' @export
#' @importFrom yaml read_yaml
#' @rdname test-helpers
#' @return the list-of-list definitions for the example `pData` returned from
#'   [example_sample_covariates()]
example_sample_covariate_definitions <- function() {
  out <- example_meta(file.path=FALSE)
  out$sample_covariate
}
