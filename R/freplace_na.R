defaults.freplace_na <- list(
  numeric = "error",
  categorical = "NA.")

#' Replaces NA's with specified values.
#'
#' Some the downstream uses of a FacilexStore can throw problems when NA's
#' are found in x or covariates, so we often want to fill in NA's with
#' non-NA markers of missing values. Note that unless specified otherwise
#' (using the `replace` and `defaults` parameters),
#'
#' Depending on the atomic type of the thing that NA's are being replaced with,
#' a default value is assumed. These can be overriden by using the `defaults`
#' parameter, or specifically by column (or list) names via the `replace`
#' parameter.
#'
#' Missing values (NA's) come up often in FacilexStores since we often use
#' them to include x from multiple xsets, which induces "ragged" (sparse)
#' covariate (px) entries. In man
#'
#' @export
#' @param x the thing that has NA's in it (a data.frame, list, or vector)
#' @param replace a named list of elements to use for custom replacement values
#' @param defaults if named elements in `x` do not appear in `replace`, you
#'   can provide default values for categories of parameters (ie.
#'   `"categorical"` or `"numeric"`), otherwise
#'   `FacileData:::defaults.freplace_na` will be used.
#' @param ignore the names of columns (or elements in a list) 
#' @return an data-replaced version of `x`
#' @examples
#' data <- data.frame(
#'   a = rnorm(10),
#'   b = letters[1:10],
#'   c = factor(LETTERS[1:10]))
#' data[3, ] <- NA
#' r1 <- freplace_na(data, list(b = "bee"), ignore = "a")
#' r2 <- freplace_na(data, list(b = "bee"), defaults = list(numeric = -Inf))
freplace_na <- function(x, replace = list(), defaults = list(),
                        ignore = character(), ...) {
  UseMethod("freplace_na", x)
}

#' @export
freplace_na.default <- function(x, replace = NULL, defaults = list(),
                                ignore = character(), ...) {
  # ignore not used here
  stopifnot(is.atomic(x))
  isna <- is.na(x)
  if (!any(isna)) return(x)

  if (is.null(defaults)) defaults <- list()
  assert_list(defaults, names = "unique")
  defaults <- c(defaults, defaults.freplace_na)
  defaults <- defaults[!duplicated(names(defaults))]

  if (is.null(replace)) {
    if (is.numeric(x)) {
      replace <- defaults[["numeric"]]
    } else {
      replace <- defaults[["categorical"]]
    }
  }

  stopifnot(is.atomic(replace) && length(replace) == 1L)

  if (is.numeric(x)) {
    if (!test_number(replace)) {
      stop("Can't replace numerics yet with anything but a number")
    }
  } else {
    replace <- as.character(replace)
    if (is.factor(x)) {
      if (!is.element(replace, levels(x))) {
        levels(x) <- c(levels(x), replace)
      }
    } else {
      x <- as.character(x)
    }
  }

  x[isna] <- replace
  x
}

#' @export
freplace_na.list <- function(x, replace = list(), defaults = list(),
                             ignore = character(), ...) {
  assert_character(ignore, null.ok = TRUE)
  for (cname in setdiff(names(x), ignore)) {
    vals <- x[[cname]]
    rep.val <- replace[[cname]]
    if (!identical(rep.val, "skip")) {
      x[[cname]] <- freplace_na(vals, replace[[cname]], defaults = defaults)
    }
  }
  x
}

#' @export
#' @method freplace_na data.frame
freplace_na.data.frame <- function(x, replace = list(), defaults = list(),
                                   ignore = character(), ...) {
  freplace_na.list(x, replace, defaults, ignore, ...)
}
