##' Creates a FacileDb connection to a test database
##'
##' TODO: This function currently uses the Atezo database as a test. We need to
##'       create a REAL testing database.
##'
##' @export
##' @param db.fn The path on the filesystem to the facile database
##' @param covdef The path on the filesystem to the yaml file that provides
##'   the covariate definitions
exampleFacileDataSet <- function(db.type=c('sqlite', 'monetdblite')) {
  fn <- system.file('extdata', 'exampleFacileDataSet', package='FacileDataSet')
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

