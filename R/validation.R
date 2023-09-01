# checkmate compliant validations functions ------------------------------------

#' Check to see if a vector is categorical (character or string)
#'
#' @export
#' @param x a vector of things
#' @param any.missing are vectors with missing values allowed? Default is `TRUE`
#' @param all.missing are vectors with missing values allowed? Default is `TRUE`
#' @param len expected length of `x`. If provided, overrides `min.len` and
#'   `max.len`. Defaults to `NULL`.
#' @param min.len minimum length for `x`
#' @param max.len maximum length for `x`
#' @param ... dots
check_categorical <- function(x, any.missing = TRUE, all.missing = TRUE,
                              len = NULL, min.len = NULL, max.len = NULL, ...) {
  e <- character()
  check_fn <- if (is.character(x)) {
    check_character
  } else if (is.factor(x)) {
    check_factor
  } else if (is.logical(x)) {
    check_logical
  } else {
    function(x) "x is not a character, factor, or logical"
  }
  check_fn(x)
}

#' @rdname check_categorical
#' @export
assert_categorical <- function(x, any.missing = TRUE, all.missing = TRUE,
                               len = NULL, min.len = NULL, max.len = NULL, ...,
                               .var.name = vname(x), add = NULL) {
  res <- check_categorical(x, any.missing = any.missing,
                           all.missing = all.missing, len = len,
                           min.len = min.len, max.len = max.len, ...)
  makeAssertion(x, res, .var.name, add)
}

#' @rdname check_categorical
#' @export
test_categorical <- function(x, ...) {
  identical(check_categorical(x, ...), TRUE)
}


#' Check if argument is a FacileDataStore
#'
#' @export
#'
#' @param x The object to check.
#' @param ... to be determined later
check_facile_data_store <- function(x, ...) {
  e <- character()
  if (!is(x, "FacileDataStore")) {
    e <- paste0("Must be of type 'FacileDataStore', not '", class(x)[1L], "'")
  }

  if (length(e)) e else TRUE
}

#' @export
#' @rdname check_facile_data_store
#' @param .var.name Name of the checked object to print in assertions. Defaults
#'   to the heuristic implemented in [checkmate::vname()].
#' @param add An [checkmate::AssertCollection()] object. Default is `NULL`.
assert_facile_data_store <- function(x, ..., .var.name = vname(x), add = NULL) {
  res <- check_facile_data_store(x, ...)
  makeAssertion(x, res, .var.name, add)
}

#' @export
#' @rdname check_facile_data_store
test_facile_data_store <- function(x, ...) {
  identical(check_facile_data_store(x, ...), TRUE)
}

#' Check if argument is a FacileDataSet
#'
#' @export
#' @inheritParams check_facile_data_store
check_facile_data_set <- function(x, ...) {
  e <- character()
  if (!is(x, "FacileDataSet")) {
    e <- paste0("Must be of type 'FacileDataSet', not '", class(x)[1L], "'")
  }
  if (!("con" %in% names(x) && is(x[["con"]], "DBIObject"))) {
    e <- c(e, "x$con is not a DBIObject")
  }
  if (!("anno.dir" %in% names(x) &&
        check_directory_exists(x[["anno.dir"]], "w"))) {
    e <- c(e, "x$anno.dir is not a valid annotation directory")
  }
  if (!("hdf5.fn" %in% names(x) &&
        check_file_exists(x[["hdf5.fn"]], "r"))) {
    e <- c(e, "x$anno.dir is not a valid HDF5 file")
  }
  if (length(e)) e else TRUE
}

#' @export
#' @rdname check_facile_data_set
assert_facile_data_set <- function(x, ..., .var.name = vname(x), add = NULL) {
  res <- check_facile_data_set(x, ..., )
  makeAssertion(x, res, .var.name, add)
}

#' @export
#' @rdname check_facile_data_set
test_facile_data_set <- function(x, ...) {
  identical(check_facile_data_set(x, ...), TRUE)
}

# checkmate-like validation functios -------------------------------------------

#' Check to see that samples are referenced correctly
#'
#' Samples have compound keys: dataset,sample_id. If we want to index into
#' them, we can either:
#'
#'   1. pass a data.frame around with dataset and sample_id columns
#'   2. pass a "loaded up" tbl_sqlite" over the sample_covariate table which
#'      has your filters of interest set
#' @export
#' @rdname assertions
assert_sample_subset <- function(x, fds = NULL, ..., .var.name = vname(x),
                                 add = NULL) {
  res <- check_sample_subset(x, fds, ...)
  makeAssertion(x, res, .var.name, add)
}

#' @export
#' @rdname assertions
check_sample_subset <- function(x, fds = NULL, ...) {
  e <- character()
  if (!(is(x, 'tbl') || is(x, 'data.frame'))) {
    e <- "Sample descriptor is not data.frame/tbl-like"
  }
  if (!has_columns(x, c('dataset', 'sample_id'), warn = FALSE)) {
    e <- c(e, "'dataset' and 'sample_id' columns required in sample descriptor")
  }
  if (length(e) == 0L && !is.null(fds)) {
    .samples <- samples(fds, .valid_sample_check = FALSE)
    bad.samples <- anti_join(x, .samples, by = c("dataset", "sample_id"),
                             copy = !same_src(.samples, x))
    bad.samples <- collect(bad.samples, n = Inf)
    nbad <- nrow(bad.samples)
    if (nbad > 0L) {
      e <- c(e, paste(nbad, "samples not found in FacileDataStore"))
    }
  }

  if (length(e)) e else TRUE
}

#' @export
#' @rdname assertions
test_sample_subset <- function(x, fds = NULL, ...) {
  identical(check_sample_subset(x, fds, ...), TRUE)
}

#' @export
#' @rdname assertions
assert_facet_descriptor <- function(x) {
  stopifnot(is_facet_descriptor(x))
  invisible(x)
}

#' @export
#' @rdname assertions
is_facet_descriptor <- function(x) {
  if (!is_sample_subset(x)) return(FALSE)
  has_columns(x, 'facet')
}

#' @section assay_feature_descriptor:
#' If .fds is provided, it must be a \code{FaclieDataSet} and these functions
#' will check to ensure that the \code{x[['assay']]} is a valid assay element
#' in \code{.fds}
#' @export
#' @rdname assertions
assert_assay_feature_descriptor <- function(x, .fds=NULL) {
  stopifnot(is_assay_feature_descriptor(x, .fds))
  invisible(x)
}

#' @export
#' @rdname assertions
is_assay_feature_descriptor <- function(x, .fds=NULL) {
  if (!(is(x, 'tbl') || is(x, 'data.frame'))) return(FALSE)
  if (!has_columns(x, c('assay', 'feature_id'))) return(FALSE)
  if (!is.null(.fds)) {
    assert_facile_data_store(.fds)
    bad.assay <- setdiff(x[['assay']], assay_names(.fds))
    if (length(bad.assay)) {
      stop("Assay(s) in assay_feature_descriptor not found: ",
           paste(bad.assay, collapse=','))
    }
  }
  TRUE
}

#' @export
#' @rdname assertions
assert_expression_result <- function(x) {
  stopifnot(is_expression_result(x))
  invisible(x)
}

#' @export
#' @rdname assertions
is_expression_result <- function(x) {
  if (!(is(x, 'tbl') || is(x, 'data.frame'))) return(FALSE)
  has_columns(x, c('dataset', 'sample_id', 'feature_id', 'value'))
}

#' @export
#' @rdname assertions
assert_sample_statistics <- function(x) {
  stopifnot(is_sample_statistics(x))
  invisible(x)
}

#' @export
#' @rdname assertions
is_sample_statistics <- function(x) {
  if (!(is(x, 'tbl') || is(x, 'data.frame'))) return(FALSE)
  has_columns(x, c('dataset', 'sample_id', 'libsize', 'normfactor'))
}

#' @export
#' @rdname assertions
assert_sample_covariates <- function(x) {
  stopifnot(is_sample_covariates(x))
  invisible(x)
}

#' @export
#' @rdname assertions
is_sample_covariates <- function(x) {
  if (!(is(x, 'tbl') || is(x, 'data.frame'))) return(FALSE)
  req.cols <- c('dataset', 'sample_id', 'variable', 'value', 'class') #, 'type')
  has_columns(x, req.cols)
}

#' @export
#' @rdname assertions
assert_columns <- function(x, req.cols) {
  stopifnot(has_columns(x, req.cols))
  invisible(x)
}

#' @export
#' @rdname assertions
has_columns <- function(x, req.cols, warn = TRUE) {
  missed <- setdiff(req.cols, colnames(x))
  any.missing <- length(missed) > 0L
  if (any.missing && warn) {
    warning("missing columns: ", paste(missed, collpase=', '), immediate.=TRUE)
  }
  !any.missing
}

#' @export
#' @rdname assertions
assert_covariate_definitions <- function(x, required = NULL) {
  stopifnot(is_covariate_definitions(x, required))
  invisible(x)
}

#' @export
#' @rdname assertions
is_covariate_definitions <- function(x, required = NULL) {
  if (!is.list(x)) return(FALSE)
  if (!is.character(names(x))) return(FALSE)
  if (!all(sapply(x, is.list))) return(FALSE)

  if (length(required)) {
    # required <- c('type', 'class', 'description', 'label')
    kosher <- sapply(x, function(y) {
      sapply(required, function(z) is.character(y[[z]]))
    }) |> t()
    if (!all(kosher)) return(FALSE)
  }

  TRUE
}

