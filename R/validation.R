##' Check to see that samples are referenced correctly
##'
##' Samples have compound keys: dataset,sample_id. If we want to index into
##' them, we can either:
##'
##'   1. pass a data.frame around with dataset and sample_id columns
##'   2. pass a "loaded up" tbl_sqlite" over the sample_covariate table which
##'      has your filters of interest set
assert_sample_subset <- function(x) {
  key.cols <- c('dataset', 'sample_id')
  stopifnot(is(x, 'tbl') || is(x, 'data.frame'))
  assert_columns(x, key.cols)

  if (is.data.frame(x)) {
    ## stopifnot(nrow(x) >= 1L)
  } else if (is(x, 'tbl_sqlite')) {
    ## stopifnot(isTRUE(as.character(x$from) == 'sample_covariate'))
    ## stopifnot(is.list(x$where) && length(x$where) >= 1L)
    ## stopifnot(is.list(x$select) && length(x$select) == 2L &&
    ##             setequal(names(x$select), c('dataset', 'sample_id')))
  }

  invisible(x)
}

assert_expression_result <- function(x) {
  stopifnot(is(x, 'tbl') || is(x, 'data.frame'))
  assert_columns(x, c('dataset', 'sample_id', 'feature_id', 'count'))
  invisible(x)
}

assert_sample_statistics <- function(x) {
  stopifnot(is(x, 'tbl') || is(x, 'data.frame'))
  assert_columns(x, c('dataset', 'sample_id', 'libsize', 'normfactor'))
  invisible(x)
}

assert_sample_covariates <- function(x) {
  stopifnot(is(x, 'tbl') || is(x, 'data.frame'))
  assert_columns(x, c('dataset', 'sample_id', 'variable', 'value'))
  invisible(x)
}

assert_columns <- function(x, req.cols) {
  missed <- setdiff(req.cols, colnames(x))
  if (length(missed)) {
    stop("missing columns: ", paste(missed, collpase=', '))
  }
  invisible(x)
}
