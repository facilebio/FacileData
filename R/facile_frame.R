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
  if (.valid_sample_check) assert_sample_subset(x, datastore)
  assert_class(datastore, "FacileDataStore")
  if (!is.null(classes)) assert_character(classes)
  class(x) <- unique(c(classes, "facile_frame", class(x)))
  set_fds(x, datastore)
}

# dplyr++ manipulation =========================================================

#' Retrieves the extra class information that might be ahead of the root
#' `"facile_frame` entry, so we can spank this back on to the outgoing result
#' @noRd
.extra_classes <- function(x, ..., .root = "facile_frame") {
  classes <- class(x)
  root.idx <- match(.root, classes)
  if (is.na(root.idx)) character() else head(classes, root.idx - 1L)
}

#' @export
#' @noRd
arrange.facile_frame <- function(.data, ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(.data), .extra_classes(.data),
                  .valid_sample_check = FALSE)
}

#' @export
#' @noRd
collect.facile_frame <- function(x, ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(x), .extra_classes(x), .valid_sample_check = FALSE)
}

#' @export
#' @noRd
distinct.facile_frame <- function(.data, ..., .keep_all = FALSE) {
  res <- NextMethod()
  as_facile_frame(res, fds(.data), .extra_classes(.data),
                  .valid_sample_check = FALSE)
}

#' @export
#' @noRd
filter.facile_frame <- function(.data, ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(.data), .extra_classes(.data),
                  .valid_sample_check = FALSE)
}

#' @export
#' @noRd
group_by.facile_frame <- function(.data, ..., add = FALSE) {
  res <- NextMethod()
  # as_facile_frame(res, fds(.data), .valid_sample_check = FALSE)
  res
}

#' @export
#' @noRd
mutate.facile_frame <- function(.data, ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(.data), .extra_classes(.data),
                  .valid_sample_check = FALSE)
}

#' @export
#' @noRd
select.facile_frame <- function(.data, ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(.data), .extra_classes(.data),
                  .valid_sample_check = FALSE)
}

#' @export
#' @noRd
subset.facile_frame <- function(x, ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(x), .extra_classes(x), .valid_sample_check = FALSE)
}

#' @export
#' @noRd
transmute.facile_frame <- function(.data, ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(.data), .extra_classes(.data),
                  .valid_sample_check = FALSE)
}

#' @export
#' @noRd
ungroup.facile_frame <- function(x, ...) {
  res <- NextMethod()
  # as_facile_frame(res, fds(x), .valid_sample_check = FALSE)
  res
}

# Joins ========================================================================

inner_join.facile_frame <- function(x, y, by = NULL, copy = FALSE,
                                    suffix = c(".x", ".y"), ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(x), .extra_classes(x), .valid_sample_check = FALSE)
}

left_join.facile_frame <- function(x, y, by = NULL, copy = FALSE,
                                   suffix = c(".x", ".y"), ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(x), .extra_classes(x), .valid_sample_check = FALSE)
}

right_join.facile_frame <- function(x, y, by = NULL, copy = FALSE,
                                    suffix = c(".x", ".y"), ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(x), .extra_classes(x), .valid_sample_check = FALSE)
}

full_join.facile_frame <- function(x, y, by = NULL, copy = FALSE,
                                   suffix = c(".x", ".y"), ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(x), .extra_classes(x), .valid_sample_check = FALSE)
}

anti_join.facile_frame <- function(x, y, by = NULL, copy = FALSE, ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(x), .extra_classes(x), .valid_sample_check = FALSE)
}

