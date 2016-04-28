##' Fetch the sample statistics for sets of samples in the warehouse
##'
##' @export
##' @param db the connection to the database
##' @param samples a data.frame or tbl_sqlite that has dataset and sample_id
##'   columns
##' @return a tbl_df or tbl_sqlite result from the sample_stats table
fetch_sample_statistics <- function(db, samples=NULL, do.collect=FALSE) {
  stopifnot(is.FacileDb(db))
  ss <- sample_stats_tbl(db)
  if (is.null(samples)) {
    return(ss)
  }
  assert_sample_subset(samples)
  internalize <- !same_src(ss, samples)
  out <- select(samples, dataset, sample_id) %>%
    distinct %>%
    inner_join(ss, ., by=c('dataset', 'sample_id'), copy=internalize,
                auto_index=internalize)
  if (do.collect) {
    out <- collect(out)
    db <- NULL
  }
  attr(out, 'db') <- db
  out
}

