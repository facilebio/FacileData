# dplyr extensions
# for facile_frames and other useful things in the tidyverse

#' Retrieves the extra class information that might be ahead of the root
#' `"facile_frame` entry, so we can spank this back on to the outgoing result
#' @noRd
.extra_classes <- function(x, ..., .root = "facile_frame") {
  classes <- class(x)
  root.idx <- match(.root, classes)
  if (is.na(root.idx)) character() else utils::head(classes, root.idx - 1L)
}

#' A version of group_map that provides a named list in return
#'
#' dplyr::group_map returns an unnamed list, but a named list is often handy
#' https://github.com/tidyverse/dplyr/issues/4223#issuecomment-469269857
#'
#' @export
#' @param .data a grouped tibble
#' @param .f a function of formula to apply to each group.
#'   If a function, it is used as is. It should have at least 2 formal arguments.
#'   If a formula, e.g. ~ head(.x), it is converted to a function.
#'   In the formula, you can use:
#'     *  `.` or `.x` to refer to the subset of rows of `.tbl` for the given
#'        group
#'     * `.y` to refer to the key, a one row tibble with one column per grouping
#'       variable that identifies the group
#' @return a list of elemnts returned by `.f` over the grouped elements in
#'   `.data`
#' @examples
#' # no names
#' iris |> group_by(Species) |> group_map(~ nrow(.x))
#' 
#' # with names
#' iris |> group_by(Species) |> group_map.(~ nrow(.x))
group_map. <- function(.data, .f, ..., .keep = FALSE) {
  lifecycle::signal_stage("experimental", "group_map.()")
  UseMethod("group_map.")
}

#' @export
#' @noRd
#' @method group_map. data.frame
group_map..data.frame <- function (
    .data, .f, ...,
    .keep = FALSE, keep = deprecated(),
    .sep = " / ") {
  assert_string(.sep)
  if (!missing(keep)) {
    lifecycle::deprecate_warn("1.0.0", "group_map(keep = )", 
                              "group_map.(.keep = )", always = TRUE)
    .keep <- keep
  }
  
  .f <- dplyr:::as_group_map_function(.f)
  chunks <- if (dplyr::is_grouped_df(.data)) {
    group_split.(.data, .keep = isTRUE(.keep), .sep = .sep)
  } else {
    group_split.(.data)
  }
  keys <- dplyr::group_keys(.data)
  group_keys <- map(seq_len(nrow(keys)), function(i) {
    keys[i, , drop = FALSE]
  })
  if (length(chunks)) {
    map2(chunks, group_keys, .f, ...)
  } else {
    structure(
      list(),
      ptype = .f(attr(chunks, "ptype"), keys[integer(0L),], ...))
  }
}

#' Split a grouped tbl, but return a named list.
#'
#' [dplyr::group_split()] splits a grouped tbl, but the splits are unnamed.
#' This spanks the names on there.
#' https://github.com/tidyverse/dplyr/issues/4223#issuecomment-469269857
#'
#' @export
#' @param .tbl a tibble to split
#' @param ... the `...` in [dplyr::group_by()]
#' @return a list of the splitted tbl, with named elements
#' @examples
#' # nonames
#' isplit <- dplyr::group_split(iris, Species)
#' # names
#' isplit <- group_split.(iris, Species)
group_split. <- function(.tbl, ..., .keep = TRUE) {
  lifecycle::signal_stage("experimental", "group_split()")
  UseMethod("group_split.")
}

#' @noRd
#' @export
#' @method group_split. data.frame
group_split..data.frame <- function(.tbl, ..., .keep = TRUE,  .sep = " / ") {
  if (dplyr::is_grouped_df(.tbl)) {
    grouped <- .tbl
  } else {
    grouped <- dplyr::group_by(.tbl, ...)
  }
  
  nms <- rlang::inject(paste(!!!dplyr::group_keys(grouped), sep = .sep))
  
  grouped |> 
    dplyr::group_split() |> 
    rlang::set_names(nms)
}

#' @export
#' @noRd
arrange.facile_frame <- function(.data, ..., .facilitate = NULL) {
  res <- NextMethod()
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(.data), "FacileDataStore"))
  }
  if (.facilitate) {
    res <- as_facile_frame(res, fds(.data), .extra_classes(.data),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
collect.facile_frame <- function(x, ..., .facilitate = NULL) {
  res <- NextMethod()
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(x), "FacileDataStore"))
  }
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
distinct.facile_frame <- function(.data, ..., .keep_all = FALSE,
                                  .facilitate = NULL) {
  res <- NextMethod()
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(.data), "FacileDataStore"))
  }
  if (.facilitate) {
    res <- as_facile_frame(res, fds(.data), .extra_classes(.data),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
filter.facile_frame <- function(.data, ..., .facilitate = NULL) {
  res <- NextMethod()
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(.data), "FacileDataStore"))
  }
  if (.facilitate) {
    res <- as_facile_frame(res, fds(.data), .extra_classes(.data),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @noRd
#' @export
#' @examples
#' xx <- samples(an_fds()) |> with_sample_covariates()
#' xx |> 
#'   group_by(cell_type) |> 
#'   summarize(avg = mean(ncells))
group_by.facile_frame <- function(.data, ..., add = FALSE,
                                  # .drop = group_by_drop_default(.data),
                                  .drop = TRUE,
                                  .facilitate = NULL) {
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
mutate.facile_frame <- function(.data, ..., .facilitate = NULL) {
  res <- NextMethod()
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(.data), "FacileDataStore"))
  }
  if (.facilitate) {
    res <- as_facile_frame(res, fds(.data), .extra_classes(.data),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
rename.facile_frame <- function(.data, ..., .facilitate = NULL) {
  res <- NextMethod()
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(.data), "FacileDataStore"))
  }
  if (.facilitate) {
    res <- as_facile_frame(res, fds(.data), .extra_classes(.data),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
select.facile_frame <- function(.data, ..., .facilitate = NULL) {
  res <- NextMethod()
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(.data), "FacileDataStore"))
  }
  if (.facilitate) {
    res <- as_facile_frame(res, fds(.data), .extra_classes(.data),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
subset.facile_frame <- function(x, ..., .facilitate = NULL) {
  res <- NextMethod()
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(.data), "FacileDataStore"))
  }
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
transmute.facile_frame <- function(.data, ..., .facilitate = NULL) {
  res <- NextMethod()
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(.data), "FacileDataStore"))
  }
  if (.facilitate) {
    res <- as_facile_frame(res, fds(.data), .extra_classes(.data),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
#' @examples
#' afds <- an_fds()
#' ssa <- samples(afds) |> with_sample_covariates()
#' ssg <- ssa |> 
#'   group_by(cell_abbrev) |> 
#'   mutate(n = n(), .before = 1L) |> 
#'   ungroup()
group_by.facile_frame <- function(
    .data,
    ...,
    .add = FALSE,
    .drop = group_by_drop_default(.data)) {
  fds. <- fds(.data)
  res <- NextMethod()
  as_facile_frame(res, fds., .extra_classes(res), .valid_sample_check = FALSE)
}

#' 
#' @export
#' @noRd
ungroup.facile_frame <- function(x, ..., .facilitate = NULL) {
  fds. <- fds(x)
  res <- NextMethod()
  as_facile_frame(res, fds., .extra_classes(res), .valid_sample_check = FALSE)
}

# Joins ========================================================================
#
# Around December 2022 (prep for dplyr 1.1?), dplyr was enforcing empty `...`
# in the *_join.data.frame functions via a call to check_dots_empty0(...).
# 
# Take a look at these references to get an idea of what I mean:
#   Issue: https://github.com/tidyverse/dplyr/issues/6599
#   PR: https://github.com/tidyverse/dplyr/pull/6605
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
upcast_ff <- function(x, downcasted, .facilitate = FALSE, ...) {
  cast_info <- attr(downcasted, "cast_info")
  assert_list(cast_info, min.len = 2L)
  assert_flag(.facilitate)
  
  original_class <- assert_string(cast_info[["original_class"]])
  down_class <- assert_string(cast_info[["down_class"]])
  if (isTRUE(class(x)[1L] == down_class) && 
      !(.facilitate && original_class == "facile_frame")) {
    class(x) <- c(original_class, class(x))
  }
  attr(x, "cast_info") <- NULL
  x
}

#' @export
#' @noRd
inner_join.facile_frame <- function(x, y, by = NULL, copy = FALSE,
                                    suffix = c(".x", ".y"), ...,
                                    .facilitate = NULL,
                                    keep = NULL) {
  # dplyr::inner_join.data.frame calls check_dots_empty0(...)
  xx <- downcast_ff(x)
  res <- inner_join(xx, y, by = by, copy = copy, suffix = suffix, ...,
                    keep = keep)
  res <- upcast_ff(res, xx)
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(x), "FacileDataStore"))
  }
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
                                   .facilitate = NULL, 
                                   keep = NULL) {
  # dplyr::left_join.data.frame calls check_dots_empty0(...)
  xx <- downcast_ff(x)
  res <- left_join(xx, y, by = by, copy = copy, suffix = suffix, ...,
                   keep = keep)
  res <- upcast_ff(res, xx)
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(x), "FacileDataStore"))
  }
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
                                    .facilitate = NULL,
                                    keep = NULL) {
  # dplyr::right_join.data.frame calls check_dots_empty0(...)
  xx <- downcast_ff(x)
  res <- right_join(xx, y, by = by, copy = copy, suffix = suffix, ...,
                    keep = keep)
  res <- upcast_ff(res, xx)
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(x), "FacileDataStore"))
  }
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
                                   .facilitate = NULL, keep = TRUE) {
  # dplyr::full_join.data.frame calls check_dots_empty0(...)
  xx <- downcast_ff(x)
  res <- full_join(xx, y, by = by, copy = copy, suffix = suffix, ...,
                   keep = keep)
  res <- upcast_ff(res, xx)
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(x), "FacileDataStore"))
  }
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
semi_join.facile_frame <- function(x, y, by = NULL, copy = FALSE, ...,
                                   .facilitate = NULL) { 
  # dplyr::semi_join.data.frame calls check_dots_empty0(...)
  xx <- downcast_ff(x)
  res <- semi_join(xx, y, by = by, copy = copy, ...)
  res <- upcast_ff(res, xx)
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(x), "FacileDataStore"))
  }
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
anti_join.facile_frame <- function(x, y, by = NULL, copy = FALSE, ...,
                                   .facilitate = NULL) { 
  # dplyr::anti_join.data.frame calls check_dots_empty0(...)
  xx <- downcast_ff(x)
  res <- anti_join(xx, y, by = by, copy = copy, ...)
  res <- upcast_ff(res, xx)
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(x), "FacileDataStore"))
  }
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

#' @export
#' @noRd
nest_join.facile_frame <- function(x, y, by = NULL, copy = FALSE, keep = NULL,
                                   name = NULL, ..., .facilitate = NULL) {
  # dplyr::nest_join.data.frame calls check_dots_empty0(...)
  xx <- downcast_ff(x)
  res <- nest_join(xx, y, by = by, copy = copy, keep = keep, name = name, ...)
  res <- upcast_ff(res, xx)
  if (is.null(.facilitate)) {
    .facilitate <- suppressWarnings(is(fds(x), "FacileDataStore"))
  }
  if (.facilitate) {
    res <- as_facile_frame(res, fds(x), .extra_classes(x),
                           .valid_sample_check = FALSE)
  }
  res
}

# bind =========================================================================
# No can do, the first param in these methods is `...`, which you can S3ize
