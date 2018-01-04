#' Fetch the sample statistics for sets of samples in the warehouse
#'
#' NOTE: this function needs the axe. It has been changed to use the
#' assay_sample_info_table, but the way we handle this with the new unhinged
#' assay needs to change.
#' @export
#' @param x A \code{FacileDataSet} object
#' @param samples a data.frame or tbl_sqlite that has dataset and sample_id
#'   columns
#' @param semi use \code{semi_join}? I've found this to be slow sometimes in
#'   SQLite for some reason
#' @param assay_name parameter added to keep old API same with new "unhinged"
#'   FacileDataSets.
#' @return a tbl_df or tbl_sqlite result from the sample_stats table
fetch_sample_statistics <- function(x, samples=NULL, semi=TRUE,
                                    assay_name='rnaseq') {
  stopifnot(is.FacileDataSet(x))
  assert_string(assay_name)
  stopifnot(assay_name %in% assay_names(x))

  # ss <- sample_stats_tbl(x)
  ss <- assay_sample_info_tbl(x) %>%
    filter(assay == assay_name) %>%
    set_fds(x)

  if (is.null(samples)) {
    out <- ss
  } else {
    ## TODO: Need to write unit tests here to exercise what we want to do with
    ##       these results when samples are provided
    samples <- assert_sample_subset(samples)
    out <- join_samples(ss, samples, semi)
  }

  set_fds(out, x)
}

