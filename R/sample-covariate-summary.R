#' Provides a tibble summary of the covariates available over a sample space.
#' 
#' This is similar to the `FacileData::summary.eav_covariates function`, but 
#' that function is downstream of a `collect()`. This function takes in a samples
#' facile_frame (**preferably a lazy_tbl/tbl_sql**) and summaries the covariates
#' availalbe there. Keeping this function in SQL space as long as possible
#' assists in the initial load and presentation of sample filters of monstrously
#' sized datasets.
#' 
#' @export
#' @param x an object (samples) to summarize covariates over
#' @param detailed When `TRUE`, summarizes the individual levels of covariates.
#'   Default is `FALSE`.
sample_covariate_summary <- function(x, ..., detailed = FALSE) {
  UseMethod("sample_covariate_summary", x)
}

#' @noRd
#' @export
#' @param db_copy This passed down to the `copy` parameter of the internal
#'   dplyr joining function (semi_join)
sample_covariate_summary.facile_frame <- function(x, ..., detailed = FALSE,
                                                  db_copy = FALSE) {
  sample_covariate_summary(
    fds(x), 
    samples = distinct(x, dataset, sample_id),
    detailed = detailed, db_copy = db_copy,
    ...)
}

#' @noRd
#' @export
sample_covariate_summary.FacileDataSet <- function(x, samples = NULL, ...,
                                                   categorical_only = TRUE,
                                                   with_levels = FALSE,
                                                   detailed = FALSE,
                                                   db_copy = FALSE,
                                                   verbose = FALSE) {
  sctbl <- sample_covariate_tbl(x)
  
  do_semi <- !is.null(samples)
  if (is.null(samples)) {
    samples <- samples(x)
  } else {
    if (!same_src(sctbl, samples)) {
      if (!db_copy) {
        stop("The samples table is not in the database. You can copy it into ",
             "the database to do a semi_join, but this can be an expensive ",
             "operation and you have to explicitly opt in to it with ",
             "`db_copy = TRUE`")
      } else {
        warning("performing a semi_join against an external data.frame that ",
                "needs to be copied can be deadly slow", immediate. = TRUE)
      }
    }
  }
  
  assert_class(samples, "facile_frame")
  assert_flag(detailed)
  assert_flag(categorical_only)
  assert_flag(with_levels)
  
  if (missing(categorical_only) && with_levels) {
    if (verbose) message("Setting categorical_only to TRUE")
    categorical_only <- TRUE
  } else if (missing(with_levels) && !categorical_only) {
    with_levels <- FALSE
  }
  
  if (categorical_only) {
    query <- filter(sctbl, .data$class == "categorical")
  } else {
    query <- sctbl
  }
  
  scols <- c("dataset", "sample_id", "variable", "value", "class", "type")
  if (do_semi) {
    query <- semi_join(query, samples, by = c("dataset", "sample_id"))
  }
  scols <- setdiff(scols, c("dataset", "sample_id"))
  if (!with_levels) scols <- setdiff(scols, "value")
  
  query <- select(query, all_of(scols))
  
  result <- distinct(query, .keep_all = TRUE)
  
  if (detailed) {
    # figure out distribution of real-valued covariates,
    # nlevels of categorical ones
    # number of samples each apply to?
  } else {
    # what to do?
  }
  
  collect(result, n = Inf)
}
