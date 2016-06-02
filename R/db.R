##' Setup the SQLite connection the TCGA database
##'
##' @export
##' @param db.path The path to the clinX database.
##' @return A \code{dplyr::src_sqlite} connection to the database.
FacileDb <- function(db.path=getOption('fatezo.dbpath', NULL),
                     cache_size=20000) {
  if (!is.character(db.path)) {
    stop("Either set options(facile.dbpath) or pass in path to sqlite.db file")
  }
  if (!file.exists(db.path)) {
    stop("Illegal path to tcga db file")
  }

  ## Update some parameters in the connection for increased speed
  ## http://stackoverflow.com/questions/1711631
  ##
  ## Explains some pragme setting:
  ##   http://www.codificar.com.br/blog/sqlite-optimization-faq/
  ##
  ## The following PRAGMA are of likely interest:
  ##   1. cache_size  https://www.sqlite.org/pragma.html#pragma_cache_size
  ##   2. page_size
  out <- src_sqlite(db.path)
  dbSendQuery(out$con, 'pragma temp_store=MEMORY;')
  dbSendQuery(out$con, sprintf('pragma cache_size=%d;', cache_size))

  out['cov.def'] <- list(NULL)
  out['cov.def'] <- getOption('fatezo.covdef')

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
