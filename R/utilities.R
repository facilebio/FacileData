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

#' Cast a data.frame into a matrix.
#'
#' This is a simple wrapper to \code{\link[reshape2]{dcast}} which converts
#' the \code{data.frame} that call would create into a matrix which uses the
#' first column as the rownames.
#'
#' @export
#' @param data a \code{data.frame} to whip into a matrix
#' @param formula the pivot formula
#' @param fun.aggregate the aggregate function
#' @param ... args passed down to \code{\link[reshape2]{dcast}}.
mcast <- function(data, formula, fun.aggregate=NULL, ...) {
  warning("Use acast, not mcast(?)")
  d <- dcast(data, formula, fun.aggregate, ...)
  set_rownames(as.matrix(d[, -1L, drop=FALSE]), d[[1L]])
}

#' bind_rows a list of data.frames that might have Surv columns
#'
#' Surv columns currently break bind_rows, so we coerce to character, bind_rows,
#' and coerce back.
#' @param df_list list of data.frames
#' @return data.frame
#' @export
bind_pdata_rows <- function(df_list) {
    ## FIXME: Deprecate in favor of dplyr::bind_rows when they fix #2457
    is_surv = lapply(df_list,
                     function(x) {
                         vapply(x, is, logical(1), "Surv")
                     })

    clean_df_list = mapply(df_list, is_surv,
                      FUN = function(x,y) {
                          x[y] = lapply(x[y], as.character)
                          x
                      }, SIMPLIFY = FALSE)

    big_df = bind_rows(clean_df_list)
    is_surv_uniq = unlist(is_surv)
    is_surv_uniq = is_surv_uniq[!duplicated(names(is_surv_uniq))]
    stopifnot(setequal(names(is_surv_uniq), colnames(big_df)))
    is_surv_unique = is_surv_uniq[colnames(big_df)]
    big_df[is_surv_unique] = lapply(big_df[is_surv_unique], eav_decode_Surv)
    big_df
}
