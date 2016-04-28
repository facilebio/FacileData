##' Setup the SQLite connection the TCGA database
##'
##' @export
##' @param db.path The path to the clinX database.
##' @return A \code{dplyr::src_sqlite} connection to the database.
FacileDb <- function(db.path=getOption('fatezo.dbpath', NULL)) {
  if (!is.character(db.path)) {
    stop("Either set options(facile.dbpath) or pass in path to sqlite.db file")
  }
  if (!file.exists(db.path)) {
    stop("Illegal path to tcga db file")
  }
  out <- src_sqlite(db.path)
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

##' Filters the samples down in a dataset to ones specified
##'
##' Tables like \code{expression} and \code{sample_covariate} house different
##' datapoints per sample, and we often want to only retreive data points over
##' a subset of samples.
##'
##' @export
##' @param x likely a \code{tbl_sqlite} object, but a \code{tbl_df}-like
##'   object should work as well.
##' @param samples a sample descriptor \code{tbl_df}-like object (likely a
##'   \code{tbl_sqlite} object) that has \code{"dataset"} and \code{"samle_id"}
##'   columns.
##' @return filtered version of \code{x} that only has the desired samples
filter_samples <- function(x, samples=NULL) {
  if (is.null(x)) {
    return(x)
  }
  assert_sample_subset(samples)
  internal <- !same_src(samples, x)

  ## I think I should be using `semi_join` here, but that is so slow I might
  ## as well be looking things up by hand
  inner_join(x, samples, by=c('dataset', 'sample_id'),
             copy=internal, auto_index=internal)
}
