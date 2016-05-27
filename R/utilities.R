##' Cast a data.frame into a matrix.
##'
##' This is a simple wrapper to \code{\link[reshape2]{dcast}} which converts
##' the \code{data.frame} that call would create into a matrix which uses the
##' first column as the rownames.
##'
##' @export
##' @param data a \code{data.frame} to whip into a matrix
##' @param formula the pivot formula
##' @param fun.aggregate the aggregate function
##' @param ... args past down to \code{\link[reshape2]{dcast}}.
mcast <- function(data, formula, fun.aggregate=NULL, ...) {
  d <- dcast(data, formula, fun.aggregate, ...)
  set_rownames(as.matrix(d[, -1L, drop=FALSE]), d[[1L]])
}

