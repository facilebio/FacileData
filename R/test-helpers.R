##' Creates a FacileDb connection to a test database
##'
##' TODO: This function currently uses the Atezo database as a test. We need to
##'       create a REAL testing database.
##'
##' @export
##' @param db.fn The path on the filesystem to the facile database
##' @param covdef The path on the filesystem to the yaml file that provides
##'   the covariate definitions
TestDb <- function(db.fn=getOption('ftest.dbpath', NULL),
                   covdef.fn=getOption('ftest.covdef')) {
  FacileDb(db.fn, covdef.fn)
}
