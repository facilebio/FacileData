####################################################################################################
### The FacileAPI: Classes from the FacileVerse must implement methods on each of these generics ###
####################################################################################################

#' The Facile API
#'
#' This is a stub for a manpage all about the FacileAPI

## General getters

#' @family FacileInterface
#' @export
organism <- function(x, ...) {
  UseMethod("organism")
}

#' @export
organism.default <- function(x, ...) {
  stop("The FacileAPI requires that a specific method be written for this type.")
}

#' @family FacileInterface
#' @export
samples <- function(x, ...) {
  UseMethod("samples")
}

#' @export
samples.default <- function(x, ...) {
  stop("The FacileAPI requires that a specific method be written for this type.")
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

#' @family FacileInterface
#' @export
fetch_sample_statistics <- function(x, samples=NULL, semi=TRUE, assay_name='rnaseq') {
  UseMethod("fetch_sample_statistics")
}

#' @export
fetch_sample_statistics.default <- function(x, samples=NULL, semi=TRUE, assay_name='rnaseq') {
  stop("The FacileAPI requires that a specific method be written for this type.")
}

#' @family FacileInterface
#' @export
assay_names <- function(x, default_first=TRUE) {
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

