#' Convenience wrapper to require specified packages
#'
#' @noRd
#' @param pkg A character vector of packages to require
#' @param quietly defaults to true
#' @param ... passed into [requireNamespace()]
reqpkg <- function(pkg, quietly = TRUE, ...) {
  assert_character(pkg)
  for (p in pkg) {
    if (!requireNamespace(p, ..., quietly = quietly)) {
      stop("'", p, "' package required, please install it.", call. = FALSE)
    }
  }
}


#' Arranges the columns of one data.frame to another
#'
#' This function is primarily used to add data to the FacileDataSet's SQLite
#' database. \code{x} is new data to add, and \code{to} is the a table of
#' the form that is expected in the database. We check that the columns of
#' \code{x} are a superset of columns in \code{x} and the matching columns are
#' all of the same class.
#'
#' @export
#' @param x a \code{data.frame} that needs to be checked and conformed
#' @param to the prototype \code{data.frame} that \code{x} needs to be aligned
#'   against.
#' @return the \code{tibble} version of \code{x} that is arranged to look
#'   like \code{to}.
conform_data_frame <- function(x, to) {
  stopifnot(is.data.frame(x))
  to <- suppressWarnings(collect(to, n=1L))
  stopifnot(is.data.frame(to))
  assert_columns(x, colnames(to))
  for (cname in colnames(to)) {
    if (!cname %in% colnames(x)) {
      stop("Expected columnt not found in target data.frame: ", cname)
    }
    p.class <- class(to[[cname]])[1L]
    x.class <- class(x[[cname]])[1L]
    if (p.class != x.class) {
      stop("Expected class `", p.class, "` for column '", cname, "', but got ",
           "`", x.class, "` instead")
    }
  }
  x <- as.tbl(x)
  x[, colnames(to)]
}

#' Set the class of an object and return the object
#'
#' @export
set_class <- function(x, .class, ...) {
  assert_character(.class)
  class(x) <- unique(c(.class, class(x)))
  x
}

#' Ensures that a vector has names for all elements if it has names for any
#'
#' If the vector is not named, it remains that way
#' @export
#' @param x an object with names
#' @return `x` with all elements either being uniquely named, or NULL
nameit <- function(x, ...) {
  if (is.null(names(x))) return(x)
  noname <- nchar(names(x)) == 0L
  names(x)[noname] <- x[noname]
  names(x) <- make.names(names(x), unique = TRUE)
  x
}
