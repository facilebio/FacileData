##' Check to see that samples are referenced correctly
##'
##' Samples have compound keys: dataset,sample_id. If we want to index into
##' them, we can either:
##'
##'   1. pass a data.frame around with dataset and sample_id columns
##'   2. pass a "loaded up" tbl_sqlite" over the sample_covariate table which
##'      has your filters of interest set
##' @export
##' @rdname assertions
assert_sample_subset <- function(x) {
  stopifnot(is_sample_subset(x))
  invisible(x)
}

##' @export
##' @rdname assertions
is_sample_subset <- function(x) {
  if (!(is(x, 'tbl') || is(x, 'data.frame'))) return(FALSE)
  if (!has_columns(x, c('dataset', 'sample_id'))) return(FALSE)
  TRUE
}

##' @export
##' @rdname assertions
assert_facet_descriptor <- function(x) {
  stopifnot(is_facet_descriptor(x))
  invisible(x)
}

##' @export
##' @rdname assertions
is_facet_descriptor <- function(x) {
  if (!is_sample_subset(x)) return(FALSE)
  has_columns(x, 'facet')
}

##' @export
##' @rdname assertions
assert_assay_feature_descriptor <- function(x) {
  stopifnot(is_assay_feature_descriptor(x))
  invisible(x)
}

##' @export
##' @rdname assertions
is_assay_feature_descriptor <- function(x) {
  if (!(is(x, 'tbl') || is(x, 'data.frame'))) return(FALSE)
  if (!has_columns(x, c('assay', 'feature_id'))) return(FALSE)
  TRUE
}

##' @export
##' @rdname assertions
assert_expression_result <- function(x) {
  stopifnot(is_expression_result(x))
  invisible(x)
}

##' @export
##' @rdname assertions
is_expression_result <- function(x) {
  if (!(is(x, 'tbl') || is(x, 'data.frame'))) return(FALSE)
  has_columns(x, c('dataset', 'sample_id', 'feature_id', 'count'))
}

##' @export
##' @rdname assertions
assert_sample_statistics <- function(x) {
  stopifnot(is_sample_statistics(x))
  invisible(x)
}

##' @export
##' @rdname assertions
is_sample_statistics <- function(x) {
  if (!(is(x, 'tbl') || is(x, 'data.frame'))) return(FALSE)
  has_columns(x, c('dataset', 'sample_id', 'libsize', 'normfactor'))
}

##' @export
##' @rdname assertions
assert_sample_covariates <- function(x) {
  stopifnot(is_sample_covariates(x))
  invisible(x)
}

##' @export
##' @rdname assertions
is_sample_covariates <- function(x) {
  if (!(is(x, 'tbl') || is(x, 'data.frame'))) return(FALSE)
  req.cols <- c('dataset', 'sample_id', 'variable', 'value', 'class', 'type')
  has_columns(x, req.cols)
}

##' @export
##' @rdname assertions
assert_columns <- function(x, req.cols) {
  stopifnot(has_columns(x, req.cols))
  invisible(x)
}

##' @export
##' @rdname assertions
has_columns <- function(x, req.cols) {
  missed <- setdiff(req.cols, colnames(x))
  if (length(missed)) {
    warning("missing columns: ", paste(missed, collpase=', '), immediate.=TRUE)
    FALSE
  } else {
    TRUE
  }
}
