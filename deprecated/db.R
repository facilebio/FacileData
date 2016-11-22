##' Setup the SQLite connection the TCGA database
##'
##' @export
##' @param db.path The path to the clinX database.
##' @return A \code{dplyr::src_sqlite} connection to the database.
FacileDb <- function(db.fn=getOption('ftest.dbpath', NULL),
                     covdef.fn=getOption('ftest.covdef'),
                     cache_size=20000) {
  if (!is.character(db.fn) || !file.exists(db.fn)) {
    stop("The path to db.fn is not legit")
  }
  if (!is.character(covdef.fn) || !file.exists(covdef.fn)) {
    stop("The path to covdef.fn is not legit")
  }

  db.type <- tail(strsplit(db.fn, '.', fixed=TRUE)[[1]], 1L)
  if (!db.type %in% c('sqlite', 'monetdblite')) {
    stop("Unrecognized database type. Must be *.sqlite or *.monetdblite")
  }

  if (db.type == 'sqlite') {
    ## Update some parameters in the connection for increased speed
    ## http://stackoverflow.com/questions/1711631
    ##
    ## Explains some pragme setting:
    ##   http://www.codificar.com.br/blog/sqlite-optimization-faq/
    ##
    ## The following PRAGMA are of likely interest:
    ##   1. cache_size  https://www.sqlite.org/pragma.html#pragma_cache_size
    ##   2. page_size
    out <- src_sqlite(db.fn)
    dbSendQuery(out$con, 'pragma temp_store=MEMORY;')
    dbSendQuery(out$con, sprintf('pragma cache_size=%d;', cache_size))
  } else if (db.type == 'monetdblite') {
    if (!require('MonetDBLite')) {
      stop("MonetDBLite required to access MonetDBLite database")
    }
    out <- src_monetdblite(db.fn)
  }

  out['cov.def'] <- list(NULL)
  out['cov.def'] <- covdef.fn

  class(out) <- c('FacileDb', class(out))
  out
}

##' @export
is.FacileDb <- function(x) {
  is(x, 'FacileDb') && 'cov.def' %in% names(x) && file.exists(x$cov.def)
}

##' @export
expression_tbl <- function(db=FacileDb()) {
  tbl(db, 'expression')
}

##' @export
sample_stats_tbl <- function(db=FacileDb()) {
  tbl(db, 'sample_stats')
}

##' @export
sample_covariate_tbl <- function(db=FacileDb()) {
  tbl(db, 'sample_covariate')
}

##' @export
gene_info_tbl <- function(db=FacileDb()) {
  tbl(db, 'gene_info')
}

##' Get/set db
##'
##' @rdname getsetdb
##' @export
##' @param x the object
##' @param db The \code{FacileDb} object
fdb <- function(x) {
 attr(x, 'db')
}

##' @rdname getsetdb
##' @export
"fdb<-" <- function(x, value) {
  UseMethod("db<-", x)
}

##' @rdname getsetdb
##' @export
"fdb<-.tbl" <- function(x, value) {
  attr(x, 'db') <- value
  x
}

##' @rdname getsetdb
##' @export
"fdb<-.data.frame" <- function(x, value) {
  attr(x, 'db') <- value
  x
}

"fdb<-.default" <- function(x, value) {
  attr(x, 'db') <- value
  x
}

##' @rdname getsetdb
##' @export
set_fdb <- function(x, value) {
  attr(x, 'db') <- value
  x
}
