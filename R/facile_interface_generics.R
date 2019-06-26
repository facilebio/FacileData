####################################################################################################
### The FacileAPI: Classes from the FacileVerse must implement methods on each of these generics ###
####################################################################################################

#' The Facile API
#'
#' This is a stub for a manpage all about the FacileAPI

#' Units of measure in an assay
#'
#' @export
#' @return string
assay_units <- function(x, assay_name, normalized = FALSE, abbreviate = FALSE,
                        ...) {
  UseMethod("assay_units", x)
}

#' @noRd
#' @export
assay_units.FacileDataStore <- function(x, assay_name = default_assay(x),
                                        normalized = FALSE, abbreviate = FALSE,
                                        ...) {
  ainfo <- assay_info(x, assay_name = assay_name)

  if (normalized) {
    out <- switch(ainfo$assay_type,
                  rnaseq = "log2(counts per million)",
                  isoseq = "log2(counts per million)",
                  normcounts = "log2(normalized counts)",
                  "values")
    if (abbreviate) {
      out <- sub("counts per million", "CPM", out)
    }
  } else {
    out <- switch(ainfo$assay_type,
                  rnaseq = "counts",
                  isoseq = "counts",
                  normcounts = "normalized counts",
                  "values")
  }

  out
}

#' Fetches assay meta information for the assays stored in a FacileDataStore
#'
#' @export
#' @param x A `FacileDataStore`
#' @param assay_name optional name of the assay to get information for
#' @return a tibble of meta information for the assays stored in `x`
assay_info <- function(x, assay_name = NULL, ...) {
  UseMethod("assay_info", x)
}

## General getters

#' @noRd
#' @family FacileInterface
#' @export
organism <- function(x, ...) {
  UseMethod("organism", x)
}

#' @export
organism.default <- function(x, ...) {
  stop("The FacileAPI requires that a specific method be written for this type.")
}

#' @family FacileInterface
#' @export
samples <- function(x, ...) {
  UseMethod("samples", x)
}

#' @export
samples.default <- function(x, ...) {
  stop("The FacileAPI requires that a specific method be written for this type.")
}

#' Get description of sample metadata columns
#'
#' Descriptions of the sample covariates can be specified in a FacileDataSet's
#' `meta.yaml` file. This function returns those.
#'
#' @export
#' @param x FacileDataTore
#' @param as.list single logical, return tibble or list
#' @return meta information about the sample covariates in `x`
covariate_definitions <- function(x, as.list = TRUE, ...) {
  UseMethod("covariate_definitions", x)
}

#' @family FacileInterface
#' @export
default_assay <- function(x, ...) {
  UseMethod("default_assay")
}

#' @export
default_assay.default <- function(x, ...) {
  stop("The FacileAPI requires that a specific method be written for this type.")
}

#' Retrieves grouping table for samples within a FacileDataSet.
#'
#' It is natural to define subgroups of samples within larger datasets.
#' This function returns grouping definitions (which we call "facets") for
#' a `FacileDataStore`.
#'
#' @family FacileInterface
#'
#' @param x An object of a class implementing the FacileInterface
#' @param name The specific facet (grouping) definition to return. Note that
#'   this parameter isn't yet used. Only one facet table was originally
#'   defined for each FacileDataSet, but we want to enable different facet
#'   definitions to be used in the future.
#' @return A `tibble` that defines the `dataset,sample_id` tuples that belong
#'   to each "facet" (group).
#' @export
facet_frame <- function(x, name = "default", ...) {
  UseMethod("facet_frame")
}

#' @family FacileInterface
#' @export
facet_frame.default <- function(x, name = "default", ...) {
  stop("The FacileAPI requires that a specific method be written for this type.")
}

## Filter

#' @family FacileInterface
#' @export
filter_features <- function(x, ...) {
  UseMethod("filter_features")
}

#' @export
filter_features.default <- function(x, ...) {
  stop("The FacileAPI requires that a specific method be written for this type.")
}

#' @family FacileInterface
#' @export
filter_samples <- function(x, ..., with_covariates = FALSE) {
  UseMethod("filter_samples")
}

#' @export
filter_samples.default <- function(x, ..., with_covariates = FALSE) {
  stop("The FacileAPI requires that a specific method be written for this type.")
}

## Fetch (<Mean Girls reference here>)

#' @family FacileInterface
#' @export
fetch_samples <- function(x, samples=NULL, assay="rnaseq", ...) {
  UseMethod("fetch_samples")
}

#' @export
fetch_samples.default <- function(x, samples = NULL, assay = "rnaseq", ...) {
  stop("The FacileAPI requires that a specific method be written for this type.")
}

#' NOTE: fetch_sample_statistics -> `fetch_assay_covariates`
#' Issue #2
#' @family FacileInterface
#' @export
fetch_sample_statistics <- function(x, samples = NULL, semi = TRUE,
                                    assay_name = NULL, ...) {
  UseMethod("fetch_sample_statistics")
}

#' Issue #2
#' @export
fetch_sample_statistics.default <- function(x, samples = NULL, semi = TRUE,
                                            assay_name = NULL) {
  stop("The FacileAPI requires that a specific method be written for this type.")
}

#' @family FacileInterface
#' @export
assay_names <- function(x, default_first = TRUE, ...) {
  UseMethod("assay_names")
}

#' @export
assay_names.default <- function(x, default_first = TRUE) {
  stop("The FacileAPI requires that a specific method be written for this type.")
}

#' @family FacileInterface
#' @export
assay_feature_info <- function(x, assay_name, feature_ids = NULL, ...) {
  UseMethod("assay_feature_info", x)
}

#' @family FacileInterface
#' @export
fetch_feature_info <- function(x, feature_type, feature_ids = NULL, ...) {
  UseMethod("fetch_feature_info", x)
}

fetch_feature_info.default <- function(x, feature_type, feature_ids = NULL, ...) {
  stop("Implement fetch_feature_info for class: ", class(x)[1L])
}

#' Append feature information columns to (feature-rows)
#'
#' @export
#' @param x a data.frame feature descriptor columns (feature_id, feature_type)
#' @return `x` fattened with the columns asked for
with_feature_info <- function(x, covariates = NULL, ...) {
  UseMethod("with_feature_info", x)
}

#' @noRd
#' @export
with_feature_info.default <- function(x, covariates = NULL, ...) {
  stop("Implement 'with_feature_info for: ", class(x)[1L])
}

#' Retrieve assay-specific information, this is where things like libsize and
#' normfactors are stored for RNA-seq data, maybe RIN score and other such
#' things?
#'
#' TODO: Need an assay_sample_info EAV table in FacileDataSet
#' @noRd
#' @family FacileInterface
#' @export
fetch_assay_covariates <- function(x, covariates = NULL,
                                   assay_name = default_assay(x), ...) {
  UseMethod("fetch_assay_covariates", x)
}

#' TODO: Need an assay_sample_info EAV table in FacileDataSet
#' @noRd
#' @family FacileInterface
#' @export
with_assay_covariates <- function(x, covariates = NULL,
                                   assay_name = default_assay(x), ...) {
  UseMethod("with_assay_covariates", x)
}


#' @family FacileInterface
#' @export
fetch_assay_data <- function(x, features, samples=NULL,
                             assay_name=default_assay(x),
                             normalized=FALSE, as.matrix=FALSE, ...,
                             subset.threshold=700, aggregate.by=NULL,
                             verbose=FALSE) {
  UseMethod("fetch_assay_data")
}

#' @export
fetch_assay_data.default <- function(x, features, samples=NULL,
                             assay_name=default_assay(x),
                             normalized=FALSE, as.matrix=FALSE, ...,
                             subset.threshold=700, aggregate.by=NULL,
                             verbose=FALSE) {
  stop("The FacileAPI requires that a specific method be written for this type.")
}

#' NOTE: is fetch_assay_score really necessary?
#' @family FacileInterface
#' @export
fetch_assay_score <- function(x, features, samples=NULL, assay_name=NULL,
                              as.matrix=FALSE, ..., subset.threshold=700) {
  UseMethod("fetch_assay_score")
}

#' @export
fetch_assay_score.default <- function(x, features, samples=NULL, assay_name=NULL,
                              as.matrix=FALSE, ..., subset.threshold=700) {

  stop("The FacileAPI requires that a specific method be written for this type.")
}

#' @family FacileInterface
#' @export
#' @rdname sample-covariates
#' @family API
fetch_sample_covariates <- function(x, samples=NULL, covariates=NULL,
                                    custom_key=Sys.getenv("USER"),
                                    with_source = FALSE, ...) {
  UseMethod("fetch_sample_covariates")
}

#' @export
#' @noRd
fetch_sample_covariates.default <- function(x, samples=NULL, covariates=NULL,
                                    custom_key=Sys.getenv("USER"), ...) {
  stop("The FacileAPI requires that a specific method be written for this type.")
}

#' @noRd
#' @rdname sample-covariates
#' @export
#' @family FacileInterface
fetch_custom_sample_covariates <- function(x, samples=NULL, covariates=NULL,
                                           custom_key=Sys.getenv("USER"),
                                           file.prefix="facile", ...) {
  UseMethod("fetch_custom_sample_covariates")
}

#' @noRd
#' @rdname sample-covariates
#' @export
fetch_custom_sample_covariates.default <- function(x, samples=NULL, covariates=NULL,
                                           custom_key=Sys.getenv("USER"),
                                           file.prefix="facile", ...) {
  stop("The FacileAPI requires that a specific method be written for this type.")
}

#' @export
#' @family FacileInterface
with_assay_data <- function(x, features, assay_name = NULL,
                            normalized = TRUE, aggregate.by = NULL,
                            spread = TRUE, with_assay_name = FALSE, ...,
                            verbose = FALSE, .fds = NULL) {
  UseMethod("with_assay_data", x)
}

with_assay_data.default <- function(x, features, assay_name = NULL,
                                    normalized = TRUE, aggregate.by = NULL,
                                    spread = TRUE, with_assay_name = FALSE, ...,
                                    verbose = FALSE, .fds = NULL) {
  stop("The FacileAPI requires a specific method be written for this type.")
}

#' Appends covariate columns to a query result
#'
#' Note that this function will force the collection of \code{x}
#'
#' @export
#' @rdname with_sample_covariates
#'
#' @importFrom stats complete.cases
#' @param x a facile sample descriptor
#' @param covariates character vector of covariate names. If \code{NULL}
#'   (default), returns all covariates, if is character and length() == 0, then
#'   this is a no-op (x is returned)
#' @param na.rm if \code{TRUE}, filters outgoing result such that only rows
#'   with nonNA values for the \code{covariates} specified here will be
#'   returned. Default: \code{FALSE}. Note that this will not check columns
#'   not specified in \code{covariates} for NA-ness.
#' @param custom_key The key to use to fetch more custom annotations over
#'   the given samples
#' @param .fds A \code{FacileDataSet} object
#' @return The facile \code{x} object, annotated with the specified covariates.
#'
#' @export
#' @family FacileInterface
#' @rdname sample-covariates
with_sample_covariates <- function(x, covariates = NULL, na.rm = FALSE,
                                   custom_key = Sys.getenv("USER"),
                                   .fds = NULL, ...) {
  UseMethod("with_sample_covariates", x)
}

#' @export
#' @family FacileInterface
with_sample_covariates.default <- function(x, covariates = NULL, na.rm = FALSE,
                                           custom_key = Sys.getenv("USER"),
                                           .fds = NULL, ...) {
  stop("The FacileAPI requires a specific method be written for this type.")
}

# Labeled ======================================================================

# Covaraites (and similar things) have can have both "labels" and "names". The
# "name" is the R-friendly name of the object/variable/etc. and the "label" is
# a human-readable version of the same, ie. "PFS" might be a name, and
# "Progression Free Survival" might be the label

#' Labeled acts like interface to reactive modules.
#'
#' Modules that implement this interface must return `label` and `name` reactive
#' elements within them.
#'
#' We use these when something (like a `assayFeatureSelect`) needs
#' a "computer friendly" name for itself (`name()`), or a more human readable
#' name (`label()`)
#'
#' @export
#' @rdname labeled
name <- function(x, ...) {
  UseMethod("name", x)
}

#' @noRd
#' @export
name.NULL <- function(x, ...) NULL

#' @noRd
#' @export
#' @rdname labeled
label <- function(x, ...) {
  UseMethod("label")
}

#' @noRd
#' @export
#' @rdname labeled
label.NULL <- function(x, ...) NULL

