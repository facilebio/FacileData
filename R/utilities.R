##' Like dcast but the first row becomes rownames of a matrix object
##' @export
mcast <- function(data, formula, fun.aggregate=NULL, ...) {
  d <- dcast(data, formula, fun.aggregate, ...)
  set_rownames(as.matrix(d[, -1L, drop=FALSE]), d[[1L]])
}

