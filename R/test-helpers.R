##' Creates a FacileDb connection to a test database
##'
##' TODO: This function currently uses the Atezo database as a test. We need to
##'       create a REAL testing database.
##'
##' @export
##' @param db.fn The path on the filesystem to the facile database
##' @param covdef The path on the filesystem to the yaml file that provides
##'   the covariate definitions
TestDb <- function(db.type=c('monetdblite', 'sqlite'),
                   datapath=getOption('ftest.datapath', NULL),
                   db.fn=getOption('ftest.dbpath', NULL),
                   covdef.fn=getOption('ftest.covdef')) {
  if (missing(db.fn)) {
    db.type <- match.arg(db.type)
    db.fn <- file.path(datapath, paste0('TcgaDb-test.', db.type))
  }

  FacileDb(db.fn, covdef.fn, cache_size=80000)
}

## This uses GenomicsTools to create a test database
createTestDb <- function(out.fn) {
  ## I am farming this out to the FacileTCGA package and creating a smaller
  ## version of the TCGA -- this is cicrular, I know.
  ##
  ## See FacileTCGA/inst/build-test-db
  stop("This is created in FacileTCGA/inst/build-test-db")
}

