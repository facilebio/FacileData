#' Regress out confounding variables from a data matrix.
#'
#' Data `x` is assumed to be log-like, and this function provides a simplified
#' interface to [limma::removeBatchEffect()].  The `batch` parameter replaces
#' `batch`, `batch2`, and `covariates`. The `design` parameter is replaced with
#' `main`. This function is mostly for use within the
#' `fetch_assay_data(..., normalized = TRUE, batch = 'something')` pipeline,
#' but refactored out here for general re-use.
#'
#' The `batch` and `main` parameters must be characters that will either
#' reference already existing columns in the `sample_info`, or be covariates
#' that can be retrieved from a FacileDataStore that is attached to the
#' sample_info facile_frame.
#'
#' We'll use these parameters to build a model.matrix with main and batch
#' effect and follow the use of `removeBatchEffect` as outlined in the post
#' linked to below to pull the design matrix apart and call the function with
#' the corresponding `design` and `covariates` parameters:
#'
#' https://support.bioconductor.org/p/83286/#83287
#'
#' Setting the `batch.scale` parameter to `TRUE` (the default), ensures that
#' the `rowMeans` of the returned data matrix are the same as the original
#' dataset.
#'
#' @export
#' @importFrom limma removeBatchEffect
#' @param x A matrix of values that needs to be corrected
#' @param sample_info a data.frame of covariate information for the data in `x`.
#'   The rows of `sample_info` are assumed to match the columns of `x`. This
#'   data.frame should have the covariates named in `batch` and `main` to use
#'   for the correction. If `sample_info` is a `facile_frame`, we will endeavor
#'   to pull any covariate named in `batch` and `main` that do not already
#'   appear in the columns of `sample_info`. Unlike limma's removeBatchEffect,
#'   we do not try to fish out the covariate values from anywhere in the
#'   "ether". They *must* be found in this data.frame.
#' @return a corrected version of the data matrix `x`.
remove_batch_effect <- function(x, sample_info, batch = NULL, main = NULL,
                                maintain.rowmeans = TRUE, ...) {
  if (is.null(batch) || !test_character(batch, min.len = 1L)) {
    warning("No batch covariates provided for correction", immediate. = TRUE)
    return(x)
  }
  assert_matrix(x, "numeric", any.missing = FALSE, min.rows = 1, min.cols = 2,
                col.names = "unique")
  assert_multi_class(sample_info, c("data.frame", "tbl"))
  if (!nrow(sample_info) == ncol(x)) {
    stop("rows in sample_info do not match columns in x")
  }

  if (is.character(main) && length(main) == 0L) main <- NULL
  if (!is.null(main)) assert_string(main)

  retrieve.covs <- setdiff(c(batch, main), colnames(sample_info))
  if (length(retrieve.covs)) {
    fds. <- fds(sample_info)
    is.facile <- test_class(fds., "FacileDataStore") &&
      test_sample_subset(sample_info)
    if (!is.facile) {
      stop("batch covariates missing in sample_info, and sample_info is not ",
           "attached to a FacileDataStore")
    }

    sample_info <- try({
      with_sample_covariates(sample_info, retrieve.covs, .fds = fds.)
    }, silent = TRUE)
    if (is(sample_info, "try-error")) {
      stop("Covariates for batch correction could not be found: ",
           paste(retrieve.covs, collapse = ","))
    }
  }

  # At this point, the rows of sample_info should match the columns of
  # `x`. Let's triple-check that, so we don't mistakenly correct batches
  # for mismatched samples.
  if (!samples_look_concordant(colnames(x), sample_info)) {
    stop("the sample_info data.frame doesn't look to match the sample ids")
  }

  # It's possible that the main and batch covariates are all singular,
  # ie. all factors with the same level, or whatever. Let's protect against
  # that before we removeBatchEffect
  batch.df <- sample_info[, c(main, batch), drop = FALSE]
  is.singular <- sapply(batch.df, function(vals) length(unique(vals)) == 1L)

  if (!is.null(main) && is.singular[main]) main <- NULL
  batch <- setdiff(batch, names(is.singular)[is.singular])

  if (length(batch)) {
    is.num <- sapply(sample_info[, batch, drop = FALSE], is.numeric)
    if (is.null(main)) {
      if (any(!is.num)) {
        cat.mats <- lapply(batch[!is.num], function(bcov) {
          # in limma we trust (code taken from limma::removeBatchEffect)
          batch. <- droplevels(as.factor(sample_info[[bcov]]))
          contrasts(batch.) <- contr.sum(levels(batch.))
          model.matrix(~ batch.)[, -1, drop = FALSE]
        })
        cat.mats <- do.call(cbind, cat.mats)
      } else {
        cat.mats <- matrix(0, nrow = nrow(sample_info), ncol = 0)
      }
      num.mats <- sample_info[, batch[is.num], drop = FALSE]
      batch.design <- cbind(cat.mats, num.mats)
      treatment.design <- matrix(1, nrow(sample_info), 1)
    } else {
      batch.formula <- paste(batch, collapse = " + ")
      des.formula <- paste("~", main, "+", batch.formula)
      des.matrix <- model.matrix(formula(des.formula), data = sample_info)
      main.cols <- c(1, grep(sprintf("^%s", main), colnames(des.matrix)))
      treatment.design <- des.matrix[, main.cols, drop = FALSE]
      batch.design <- des.matrix[, -(main.cols), drop = FALSE]
    }
    if (maintain.rowmeans) {
      batch.design <- scale(batch.design)
    }
    x <- removeBatchEffect(x, design = treatment.design,
                           covariates = batch.design)
  }
  x
}
