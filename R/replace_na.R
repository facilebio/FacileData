defaults.freplace_na <- list(
  numeric = "error",
  categorical = "NA.")

#' Replaces NA's with specified values.
#'
#' Some the downstream uses of a FacileDataStore can throw problems when NA's
#' are found in data or covariates, so we often want to fill in NA's with
#' non-NA markers of missing values. Note that unless specified otherwise
#' (using the `replace` and `defaults` parameters),
#'
#' Depending on the atomic type of the thing that NA's are being replaced with,
#' a default value is assumed. These can be overriden by using the `defaults`
#' parameter, or specifically by column (or list) names via the `replace`
#' parameter.
#'
#' Missing values (NA's) come up often in FacileDataStores since we often use
#' them to include data from multiple datasets, which induces "ragged" (sparse)
#' covariate (pData) entries. In man
#'
#' @export
#' @param data the thing that has NA's in it (a data.frame or vector)
#' @param replace a named list of elements to use for custom replacement values
#' @param defaults if named elements in `data` do not appear in `replace`, you
#'   can provide default values for categories of parameters (ie.
#'   `"categorical"` or `"numeric"`), otherwise
#'   FacileData:::defaults.freplace_na will be used.
#' @return an NA-replaced version of `data`
#' @examples
#' data <- data.frame(
#'   a = rnorm(10),
#'   b = letters[1:10],
#'   c = factor(LETTERS[1:10]))
#' data[3, ] <- NA
#' r1 <- freplace_na(data, list(b = "bee"), ignore = "a")
#' r2 <- freplace_na(data, list(b = "bee"), defaults = list(numeric = -Inf))
freplace_na <- function(data, replace = list(), defaults = list(),
                        ignore = character(), ...) {
  UseMethod("freplace_na")
}

#' @export
freplace_na.default <- function(data, replace = NULL, defaults = list(),
                                ignore = character(), ...) {
  # ignore not used here
  stopifnot(is.atomic(data))
  isna <- is.na(data)
  if (!any(isna)) return(data)

  if (is.null(defaults)) defaults <- list()
  assert_list(defaults, names = "unique")
  defaults <- c(defaults, defaults.freplace_na)
  defaults <- defaults[!duplicated(names(defaults))]

  if (is.null(replace)) {
    if (is.numeric(data)) {
      replace <- defaults[["numeric"]]
    } else {
      replace <- defaults[["categorical"]]
    }
  }

  stopifnot(is.atomic(replace) && length(replace) == 1L)

  if (is.numeric(data)) {
    if (!test_number(replace)) {
      stop("Can't replace numerics yet with anything but a number")
    }
  } else {
    replace <- as.character(replace)
    if (is.factor(data)) {
      if (!is.element(replace, levels(data))) {
        levels(data) <- c(levels(data), replace)
      }
    } else {
      data <- as.character(data)
    }
  }

  data[isna] <- replace
  data
}

#' @export
#' @method freplace_na data.frame
freplace_na.data.frame <- function(data, replace = list(), defaults = list(),
                                   ignore = character(), ...) {
  assert_character(ignore, null.ok = TRUE)
  for (cname in setdiff(colnames(data), ignore)) {
    vals <- data[[cname]]
    rep.val <- replace[[cname]]
    if (!identical(rep.val, "skip")) {
      data[[cname]] <- freplace_na(vals, replace[[cname]], defaults = defaults)
    }
  }
  data
}
