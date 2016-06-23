##' Fetch the sample statistics for sets of samples in the warehouse
##'
##' @export
##' @param db the connection to the database
##' @param samples a data.frame or tbl_sqlite that has dataset and sample_id
##'   columns
##' @return a tbl_df or tbl_sqlite result from the sample_stats table
fetch_sample_statistics <- function(db, samples=NULL, semi=TRUE) {
  stopifnot(is.FacileDb(db))
  ss <- sample_stats_tbl(db)

  if (is.null(samples)) {
    out <- ss
  } else {
    ## TODO: Need to write unit tests here to exercise what we want to do with
    ##       these results when samples are provided
    samples <- assert_sample_subset(samples)
    out <- join_samples(ss, samples, semi)
  }

  set_fdb(out, db)
}

