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
samples.facile_frame <- function(x, ...) {
  reqd <- setdiff(c("dataset", "sample_id"), colnames(x))
  if (length(reqd)) {
    stop("Missing required columns: ", paste(reqd, sep = ","))
  }
  distinct(x, dataset, sample_id, .keep_all = TRUE)
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

# dplyr++ manipulation =========================================================

#' Retrieves the extra class information that might be ahead of the root
#' `"facile_frame` entry, so we can spank this back on to the outgoing result
#' @noRd
.extra_classes <- function(x, ..., .root = "facile_frame") {
  classes <- class(x)
  root.idx <- match(.root, classes)
  if (is.na(root.idx)) character() else utils::head(classes, root.idx - 1L)
}

#' @export
#' @noRd
arrange.facile_frame <- function(.data, ..., .facilitate = TRUE) {
  res <- NextMethod()
  if (.facilitate) {
    res <- as_facile_frame(res, fds(.data), .extra_classes(.data),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
collect.facile_frame <- function(x, ..., .facilitate = TRUE) {
  res <- NextMethod()
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
distinct.facile_frame <- function(.data, ..., .keep_all = FALSE,
                                  .facilitate = TRUE) {
  res <- NextMethod()
  if (.facilitate) {
    res <- as_facile_frame(res, fds(.data), .extra_classes(.data),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
filter.facile_frame <- function(.data, ..., .facilitate = TRUE) {
  res <- NextMethod()
  if (.facilitate) {
    res <- as_facile_frame(res, fds(.data), .extra_classes(.data),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
group_by.facile_frame <- function(.data, ..., add = FALSE,
                                  # .drop = group_by_drop_default(.data),
                                  .drop = TRUE,
                                  .facilitate = TRUE) {
  # fds. <- fds(.data)
  # groups <- group_by_prepare(.data, ..., add = add)
  # groups[["data"]] <- lapply(groups[["data"]], set_fds, fds.)
  # # set_fds(grouped_df(groups$data, groups$group_names, .drop), fds.)
  # set_fds(grouped_ff(groups$data, groups$group_names, .drop), fds.)

  # if (.warn) {
  #   warning("group_by with facile_frames may get weird")
  # }
  res <- NextMethod()
  # set_fds(res, fds.)
  res
}

# #' @export
# #' @noRd
# ungroup.

#' @export
#' @noRd
mutate.facile_frame <- function(.data, ..., .facilitate = TRUE) {
  res <- NextMethod()
  if (.facilitate) {
    res <- as_facile_frame(res, fds(.data), .extra_classes(.data),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
select.facile_frame <- function(.data, ..., .facilitate = TRUE) {
  res <- NextMethod()
  if (.facilitate) {
    res <- as_facile_frame(res, fds(.data), .extra_classes(.data),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
subset.facile_frame <- function(x, ..., .facilitate = TRUE) {
  res <- NextMethod()
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
transmute.facile_frame <- function(.data, ..., .facilitate = TRUE) {
  res <- NextMethod()
  if (.facilitate) {
    res <- as_facile_frame(res, fds(.data), .extra_classes(.data),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
ungroup.facile_frame <- function(x, ..., .facilitate = TRUE) {
  res <- NextMethod()
  # as_facile_frame(res, fds(x), .valid_sample_check = FALSE)
  res
}

# Joins ========================================================================
#
# Around December 2022 (prep for dplyr 1.1?), dplyr was enforcing empty `...`
# in the *_join.data.frame functions via a call to check_dots_empty0(...).
# 
# This complicates the implementation of the .facile_frame joins because we
# can't simply call `NextMethod()` since our functions accept a .facilitate
# parameter, which gets passed down to the "next" function in the stack, and
# eventually it will fail the check_dots_empty0() check in
# the dplyr::*_join.data.frame() funcitons.
#
# To get around this we temporarily downcast the facile_frame `x` and use some
# heuristics to try and upcast again after the fact. I fear this is going to
# cause all sorts of problems, but let's see ...

#' @noRd
downcast_ff <- function(x, ...) {
  oclass <- class(x)[1L]
  class(x) <- class(x)[-1]
  attr(x, "cast_info") <- list(
    original_class = oclass,
    down_class = class(x)[1L])
  x
}

#' @noRd
upcast_ff <- function(x, downcasted, ...) {
  cast_info <- attr(downcasted, "cast_info")
  assert_list(cast_info, min.len = 2L)
  original_class <- assert_string(cast_info[["original_class"]])
  down_class <- assert_string(cast_info[["down_class"]])
  if (isTRUE(class(x)[1L] == down_class)) {
    class(x) <- c(original_class, class(x))
  }
  x
}

#' @export
#' @noRd
inner_join.facile_frame <- function(x, y, by = NULL, copy = FALSE,
                                    suffix = c(".x", ".y"), ...,
                                    .facilitate = TRUE,
                                    keep = NULL) {
  # dplyr::inner_join.data.frame calls check_dots_empty0(...)
  xx <- downcast_ff(x)
  res <- inner_join(xx, y, by = by, copy = copy, suffix = suffix, ...,
                    keep = keep)
  res <- upcast_ff(res, xx)
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
left_join.facile_frame <- function(x, y, by = NULL, copy = FALSE,
                                   suffix = c(".x", ".y"), ...,
                                   .facilitate = TRUE, 
                                   keep = NULL) {
  # dplyr::left_join.data.frame calls check_dots_empty0(...)
  xx <- downcast_ff(x)
  res <- left_join(xx, y, by = by, copy = copy, suffix = suffix, ...,
                   keep = keep)
  res <- upcast_ff(res, xx)
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
right_join.facile_frame <- function(x, y, by = NULL, copy = FALSE,
                                    suffix = c(".x", ".y"), ...,
                                    .facilitate = TRUE,
                                    keep = NULL) {
  # dplyr::right_join.data.frame calls check_dots_empty0(...)
  xx <- downcast_ff(x)
  res <- right_join(xx, y, by = by, copy = copy, suffix = suffix, ...,
                   keep = keep)
  res <- upcast_ff(res, xx)
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
full_join.facile_frame <- function(x, y, by = NULL, copy = FALSE,
                                   suffix = c(".x", ".y"), ...,
                                   .facilitate = TRUE, keep = TRUE) {
  # dplyr::full_join.data.frame calls check_dots_empty0(...)
  xx <- downcast_ff(x)
  res <- full_join(xx, y, by = by, copy = copy, suffix = suffix, ...,
                    keep = keep)
  res <- upcast_ff(res, xx)
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
semi_join.facile_frame <- function(x, y, by = NULL, copy = FALSE, ...,
                                   .facilitate = TRUE) { 
  # dplyr::semi_join.data.frame calls check_dots_empty0(...)
  xx <- downcast_ff(x)
  res <- semi_join(xx, y, by = by, copy = copy, ...)
  res <- upcast_ff(res, xx)
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
anti_join.facile_frame <- function(x, y, by = NULL, copy = FALSE, ...,
                                   .facilitate = TRUE) { 
  # dplyr::anti_join.data.frame calls check_dots_empty0(...)
  xx <- downcast_ff(x)
  res <- anti_join(xx, y, by = by, copy = copy, ...)
  res <- upcast_ff(res, xx)
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
nest_join.facile_frame <- function(x, y, by = NULL, copy = FALSE, keep = NULL,
                                   name = NULL, ..., .facilitate = TRUE) {
  # dplyr::nest_join.data.frame calls check_dots_empty0(...)
  xx <- downcast_ff(x)
  res <- nest_join(xx, y, by = by, copy = copy, keep = keep, name = name, ...)
  res <- upcast_ff(res, xx)
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

# bind =========================================================================
# No can do, the first param in these methods is `...`, which you can S3ize
