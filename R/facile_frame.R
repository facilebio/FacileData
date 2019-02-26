#' Converts a normal tibble/data.frame to a facile_frame
#'
#' @param x a sample-like descriptor
#' @param .fds the FacileDataStore tied to x
#' @param ... dots
#' @param .valid_sample_check If `TRUE` (default), will check if `x` is a valid
#'   subset of the FacileDataStore `.fds`. Internal functions may set this to
#'   `TRUE` to avoid the check to (1) save time; and (2) save infinite
#'   recursion in the call to `assert_sample_subset`.
as_facile_frame <- function(x, .fds = fds(x), ..., .valid_sample_check = TRUE) {
  if (.valid_sample_check) assert_sample_subset(x, .fds)
  assert_class(.fds, "FacileDataStore")
  class(x) <- unique(c("facile_frame", class(x)))
  set_fds(x, .fds)
}

# dplyr++ manipulation =========================================================

#' @export
#' @noRd
collect.facile_frame <- function(x, ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(x), .valid_sample_check = FALSE)
}

#' @export
#' @noRd
distinct.facile_frame <- function(.data, ..., .keep_all = FALSE) {
  res <- NextMethod()
  as_facile_frame(res, fds(.data), .valid_sample_check = FALSE)
}

#' @export
#' @noRd
filter.facile_frame <- function(.data, ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(.data), .valid_sample_check = FALSE)
}

#' @export
#' @noRd
select.facile_frame <- function(.data, ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(.data), .valid_sample_check = FALSE)
}

#' @export
#' @noRd
subset.facile_frame <- function(x, ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(.data), .valid_sample_check = FALSE)
}

# Joins ========================================================================

inner_join.facile_frame <- function(x, y, by = NULL, copy = FALSE,
                                    suffix = c(".x", ".y"), ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(x), .valid_sample_check = FALSE)
}

left_join.facile_frame <- function(x, y, by = NULL, copy = FALSE,
                                   suffix = c(".x", ".y"), ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(x), .valid_sample_check = FALSE)
}

right_join.facile_frame <- function(x, y, by = NULL, copy = FALSE,
                                    suffix = c(".x", ".y"), ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(x), .valid_sample_check = FALSE)
}

full_join.facile_frame <- function(x, y, by = NULL, copy = FALSE,
                                   suffix = c(".x", ".y"), ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(x), .valid_sample_check = FALSE)
}

anti_join.facile_frame <- function(x, y, by = NULL, copy = FALSE, ...) {
  res <- NextMethod()
  as_facile_frame(res, fds(x), .valid_sample_check = FALSE)
}

