#' Internal helper functions to normalize assay data into log2 space.
#'
#' This only works for the supporte assay_types. These functions should not
#' be exported.
#'
#' @param x A matrix of raw/unnormalized assay data retrieved from
#'   within the `fetch_assay_data()` itself.
#' @param features a feature descriptor data.frame that includes the
#'   feature_id's of the rows in `x`, as well as the assay name/type they
#'   were pulled from. We assert that all features come from the same assay
#'   type, and the rows here match 1:1 the rows in `x`.
#' @param samples a sample descriptor for the columns in `x`. Rows here
#'   should match columns in `x` 1:1.
#' @param batch,main paramters sent to [remove_batch_effects()] after
normalize_assay_data <- function(x, features, samples, batch = NULL,
                                 log = TRUE, prior.count = 0.1,
                                 main = NULL, verbose = FALSE, ...) {
  stopifnot(
    nrow(x) == nrow(features),                 # x and features are concordant
    all(rownames(x) == features$feature_id),   # x and features are concordant
    ncol(x) == nrow(samples),                  # x and samples are concordant
    all(colnames(x) == samples$samid),         # x and samples are concordant
    is.character(features$assay_type),         # only working with 1 assay_type
    length(unique(features$assay_type)) == 1L) # only working with 1 assay_type

  # Retrieves the log2-like assay-normalized version of the assay data
  atype <- features$assay_type[1L]
  norm.fn.name <- paste0("normalize_assay_matrix.", atype)
  normfn <- getFunction(norm.fn.name)
  out <- normfn(x, features, samples, log = log, prior.count = prior.count, ...)

  # Now batch correct data if desired (and if explicitly log2-like)
  if (test_character(batch, min.len = 1L)) {
    if (!log) stop("Do not know how to batch correct un-logged data")
    out <- remove_batch_effect(out, samples, batch, main, ...)
  }

  out
}

# Assay-specific normalization methods -----------------------------------------

#' @noRd
#' @importFrom edgeR cpm
normalize_assay_matrix.rnaseq <- function(x, features, samples,
                                          log = TRUE, prior.count = 0.1,
                                          verbose = FALSE, ...) {
  # libsize and normfactor are currently extracted from the database
  # but this will change when we suport sample_assay covariates
  assert_numeric(samples[["libsize"]])
  assert_numeric(samples[["normfactor"]])

  # the default libsize and normfactor can be overridden by the user if they
  # passed in a samples descriptor with the DGEList$samples `lib.size` and
  # `norm.factors` columns appended to it. If not, then we use the default
  # values to calculate the effective library size for the samples
  if (test_numeric(samples[["lib.size"]])) {
    lsize <- samples[["lib.size"]]
  } else  {
    lsize <- samples[["libsize"]]
  }
  if (test_numeric(samples[["norm.factors"]])) {
    nf <- samples[["norm.factors"]]
  } else {
    nf <- samples[["normfactor"]]
  }
  edgeR::cpm(x, lsize * nf, log = log, prior.count = prior.count)
}

#' @noRd
normalize_assay_matrix.isoseq <- normalize_assay_matrix.rnaseq

#' someone processed their data with salmon or kallisto and wanted to store
#' tpm. Normalizing this is just log2(val + prior.count)
#' @noRd
normalize_assay_matrix.tpm <- function(x, features, samples,
                                       log = TRUE, prior.count = 0.1,
                                       verbose = FALSE, ...) {
  out <- x + prior.count
  if (log) {
    out <- log2(out)
  }
  out
}

#' `lognorm` assay data already normalized, nothing more to do
#' @noRd
normalize_assay_matrix.lognorm <- function(x, features, samples, ...) {
  x
}

#' no internal normalization for `qpcrct` assay data yet, needs to be
#' normalized externally and saved her -- essentialyl like `lognorm` data
#' @noRd
normalize_assay_matrix.qpcrct <- function(x, features, samples, ...) {
  x
}

#' no internal normalization for `qpcrdct` assay data yet, needs to be
#' normalized externally and saved her -- essentialyl like `lognorm` data
#' @noRd
normalize_assay_matrix.qpcrdct <- function(x, features, samples, ...) {
  x
}

