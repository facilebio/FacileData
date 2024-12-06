#' Converts a normal tibble/data.frame to a facile_frame
#'
#' @export
#'
#' @param x a sample-like descriptor
#' @param datastore the FacileDataStore tied to x
#' @param classes more classes to append to the outgoing object. The
#'   `"facile_frame"` class entry is always the last one of the bunch.
#' @param ... dots
#' @param .valid_sample_check If `TRUE` (default), will check if `x` is a valid
#'   subset of the FacileDataStore `.fds`. Internal functions may set this to
#'   `TRUE` to avoid the check to (1) save time; and (2) save infinite
#'   recursion in the call to `assert_sample_subset`.
as_facile_frame <- function(x, datastore = fds(x), classes = NULL, ...,
                            .valid_sample_check = TRUE) {
  assert_class(datastore, "FacileDataStore")
  assert_multi_class(x, c("tbl", "data.frame"))

  .valid_sample_check <- .valid_sample_check &&
    has_columns(x, c("dataset", "sample_id"), warn = FALSE)
  if (.valid_sample_check) {
    assert_sample_subset(x, datastore)
  }

  if (!is.null(classes)) assert_character(classes)
  if (!is(x, "facile_frame")) {
    classes <- c(classes, "facile_frame")
  }
  classes <- setdiff(classes, class(x))
  if (length(classes)) {
    class(x) <- unique(c(classes, class(x)))
  }
  set_fds(x, datastore)
}

#' @noRd
#' @export
samples.facile_frame <- function(x, ..., dropped = FALSE, .keep_all = TRUE) {
  reqd <- setdiff(c("dataset", "sample_id"), colnames(x))
  if (length(reqd)) {
    stop("Missing required columns: ", paste(reqd, sep = ","))
  }
  if (dropped) {
    out <- attr(x, "samples_dropped")
    if (is.null(out)) {
      out <- dplyr::tibble(dataset = character(), sample_id = character())
    } else {
      out <- distinct(out, dataset, sample_id, .keep_all = .keep_all)
    }
    out <- set_fds(out, fds(x))
  } else {
    out <- distinct(x, dataset, sample_id, .keep_all = .keep_all)
  }
  out
}

#' @noRd
#' @export
organism.facile_frame <- function(x, ...) {
  organism(fds(x), ...)
}

#' For some reason, as of RStudio 1.3 preview, facile_frames became very slow
#' to render inline in an Rmd document when the FacileDataStore came along
#' with it.
#'
#' TODO: We may want to add the possible covariates you can pull from the
#' database in the printed output, just like tibbles show a list of column
#' names that are not shown in the (horizontal) output.
#'
#' @noRd
#' @export
print.facile_frame <- function(x, ..., n = NULL, width = NULL, n_extra = NULL) {
  x <- set_fds(x, NULL)
  class(x) <- setdiff(class(x), "facile_frame")
  NextMethod()
}
