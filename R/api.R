# Defines generic used across the facileverse.
#
# This is split into two parts:
# 1. Generic functions for a wide variety of objects (like `facilitate`)
# 2. The FacileData API itself, for querying and retrieving data from
#    multi-assay genomics data stores.

# Universal Funcitons ----------------------------------------------------------

#' Materialize a Bioconductor assay container from some facile object.
#'
#' Most often, this will be from some facile_frame to create a Bioconductor
#' assay container object, but this function can be overloaded for other
#'  purposes.
#'
#' The FacileAnalysis package, for example, uses this function to materialize
#' bioconductor objects of different flavors from different analysis results,
#' ie. a DGEList, or perhaps a limma fit object, etc.
#'
#' @export
#' @param x A facile object
biocbox <- function(x, ...) {
  UseMethod("biocbox", x)
}

#' Converts an arbitrary object into one that works in the facile ecosystem.
#'
#' There will be many times when the particular analysis you want to conduct
#' is not well supported in the facileverse. In this case, we will endeavor
#' to implement ways for you to take these results and bring them back into
#' the facile ecosystem so that you can benefit from the interactivity provided
#' therein.
#'
#' We'll want to define `facilitate()` over a wide variety of objects. For
#' instance:
#'
#' * `facilitate(a_DGElist)` would convert an [edgeR::DGEList()] object into
#'   a `FacileDGEList`, which is just the same DGEList that implements the
#'   FacileData API. This is a work in progress and will be implemented in the
#'   FacileBioc package.
#'
#' * You might perform a differential expression analysis using standard a
#'   standard limma pipeline, but you'll want to be able to drop this result
#'   into the facile ecosystem provided in the FacileAnalysis package.
#'   The particulars of this `faciltate()` implementation would be defined in
#'   the FacileAnalysis package, and migth look something like this:
#'
#'    ```
#'    fit <- eBayes(lmFit(elist, design))
#'    limma.res <- topTable(fit, coef = "something", n = Inf)
#'    facile.res <- facilitate(elist, fit, limma.res)
#'    ```
#'
#' It's not clear how well well we'll be able to do this, or if this is even
#' the right way to do it, but we'll need to do something.
#'
#' @seealso https://github.com/facilebio/FacileBiocData
#'
#' @export
#' @param x A non-facile object that we want to bring into the facile ecosystem
#' @param ... we're going to need a lot of flexibility in the implementation of
#'   this function for different types of analyses
#' @return A facile-subclass of `x` that can take advantage of the interactive
#'   facile ecosystem.
facilitate <- function(x, ...) {
  UseMethod("facilitate", x)
}

# fds() getters/setters ........................................................

#' Get or set the FacileDataStore for an object
#'
#' FacileDataStores are passed along with most every object generated from
#' functions in the facilebio universe. This makes it convenient to dig back
#' into a large genomics objects to retrieve data from "slim" results, like
#' a sample covariate data.frame.
#'
#' @rdname fds
#' @export
#' @param x the object
#' @param value The \code{FacileDataStore} object
fds <- function(x, ...) {
  UseMethod("fds", x)
}

#' @export
#' @rdname fds
fds.FacileDataStore <- function(x) {
  return(x)
}


#' @export
#' @rdname fds
fds.default <- function(x, ...) {
  out <- attr(x, 'fds')
  if (is.null(out)) {
    warning("No FacileDataStore found in x (", class(x)[1L], ")",
            immediate.=TRUE)
  }
  out
}

#' @rdname fds
#' @export
"fds<-" <- function(x, value) {
  UseMethod("fds<-", x)
}

#' @rdname fds
#' @export
"fds<-.tbl" <- function(x, value) {
  attr(x, 'fds') <- value
  x
}

#' @rdname fds
#' @export
"fds<-.data.frame" <- function(x, value) {
  attr(x, 'fds') <- value
  x
}

#' @export
#' @noRd
"fds<-.default" <- function(x, value) {
  attr(x, 'fds') <- value
  x
}

#' @rdname fds
#' @export
set_fds <- function(x, value) {
  attr(x, "fds") <- value
  x
}

# Labeled ----------------------------------------------------------------------

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

# FacileData API ---------------------------------------------------------------

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
#' @return a tibble of meta information for the assays stored in `x`, with these
#'   columns:
#'
#'   * `assay <chr>`: Name of the assay
#'   * `assay_type <chr>`: `"rnaseq"`, `"lognorm"`, etc. Look at
#'     `FacileData:::.assay.types` vector
#'   * `feature_type <chr>`: A string from `FacileData:::.feature.types`, ie.
#'      `"ensgid"`, `"entrez"`, `"custom"`, etc.
#'   * `description <chr>`: string description
#'   * `nfeatures <int>`: number of features we have info for
#'   * `storage_mode <chr>`: `"integer"`, `"numeric"`
assay_info <- function(x, assay_name = NULL, ...) {
  UseMethod("assay_info", x)
}

#' Utility functions to get row and column indices of rnaseq hdf5 files.
#'
#' This is called to get things like hdf5_index and scaling factors for
#' the samples in a given assay.
#'
#' @export
#' @param x \code{FacileDataStore}
#' @param assay_name the name of the assay
#' @param samples a sample descriptor
#' @return an updated version of \code{samples} decorated with hd5_index,
#'   scaling factors, etc. Note that rows in \code{samples} that do not appear
#'   in \code{assay_name} will be returnd here with NA values for hd5_index and
#'   such.
assay_sample_info <- function(x, assay_name, samples = NULL, ...) {
  UseMethod("assay_sample_info", x)
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
  assay_names(x, ...)[1L]
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

#' Returns a table of information about the features (from an assay, or ...)
#'
#' @export
#' @param x a facile object
#' @return a tibble with containing feature_id, feature_type, and whatever other
#'   columns are appropriate given `x`
features <- function(x, ...) {
  UseMethod("features", x)
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

#' @noRd
#' @export
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


#' Fetch assay data from single assay of choice
#'
#' The `(fetch|with)_assay_data` functions are some of the main workhose
#' functions of the facile ecosystem. These calls enable you to retrieve
#' raw and noramlized assay data from a FacileData container.
#'
#' `fetch_assay_data(x, ...)` will return the data in long form.
#' `with_assay_data(x, ...)` is most typically used when you already have
#' a dataset `x` (a `facile_frame`) that you want to decorate with more assay
#' data. The assay data asked for will be appended on to `x` in wide format.
#' Because `fetch` is (most often) used at a lower level of granularity,
#' `normalize` is by default set to `FALSE`, while it is set to `TRUE` in
#' `with_assay_data`.
#'
#' @section Removing Batch Effects:
#' When normalized data is returned, we assume these data are log-like, and you
#' have the option to regress out batch effects using our
#' [remove_batch_effect()] wrapper to [limma::removeBatchEffect()].
#'
#' @rdname fetch_assay_data
#' @export
#'
#' @inheritParams remove_batch_effect
#'
#' @param x A `FacileDataSrote` object, or `facile_frame`
#' @param features a feature descriptor (data.frame with assay and feature_id
#'   columms)
#' @param samples a sample descriptor to specify which samples to return data
#'   from.
#' @param assay_name the name of the assay to fetch data from. Defaults to the
#'   value of [default_assay()] for `x`. Must be a subset of `assay_names(x)`.
#' @param normalized return normalize or raw data values, defaults to `FALSE`.
#'   This is only really "functional" for for `assay_type = "rnaseq"` types
#'   of assays, where the normalized data is log2(CPM). These values can
#'   be tweaked with `log = (TRUE|FALSE)` and `prior.count` parameters, which
#'   can passed down internally to (eventually) [edgeR::cpm()].
#' @param as.matrix by default, the data is returned in a long-form tbl-like
#'   result. If set to `TRUE`, the data is returned as a matrix.
#' @param ... parameters to pass to normalization methods
#' @param subset.threshold sometimes fetching all the genes is faster than
#'   trying to subset. We have to figure out why that is, but I've previously
#'   tested random features of different lengths, and around 700 features was
#'   the elbow.
#' @param aggregate.by do you want individual level results or geneset
#'   scores? Use 'ewm' for eigenWeightedMean, and that's all.
#' @return A `tibble` (lazy or not) with assay data.
#' @examples
#' samples <- exampleFacileDataSet() %>%
#'   filter_samples(indication == "BLCA", sample_type == "tumor")
#' features <- c(PRF1='5551', GZMA='3001', CD274='29126')
#' dat <- with_assay_data(samples, features, normalized = TRUE, batch = "sex")
#' dat <- with_assay_data(samples, features, normalized = TRUE,
#'                        batch = c("sex", "stage"))
#' dat <- with_assay_data(samples, features, normealized = TRUE,
#'                        batch = c("sex", "stage"), main = "sample_type")
fetch_assay_data <- function(x, features, samples = NULL,
                             assay_name = ndefault_assay(x),
                             normalized = FALSE, batch = NULL, main = NULL,
                             as.matrix = FALSE, ...,
                             subset.threshold=700, aggregate = FALSE,
                             aggregate.by = "ewm", verbose=FALSE) {
  UseMethod("fetch_assay_data")
}

#' @export
fetch_assay_data.default <- function(x, features, samples=NULL,
                             assay_name=default_assay(x),
                             normalized=FALSE, batch = NULL, main = NULL,
                             as.matrix=FALSE, ...,
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
                            normalized = TRUE, aggregate = FALSE,
                            aggregate.by = "ewm", spread = TRUE,
                            with_assay_name = FALSE, ..., verbose = FALSE,
                            .fds = NULL) {
  UseMethod("with_assay_data", x)
}

#' @noRd
#' @export
with_assay_data.default <- function(x, features, assay_name = NULL,
                                    normalized = TRUE, aggregate = FALSE,
                                    aggregate.by = "ewm",
                                    spread = TRUE, with_assay_name = FALSE, ...,
                                    verbose = FALSE, .fds = NULL) {
  stop("The FacileAPI requires a specific method be written for this type.")
}

#' Appends covariate columns to a query result
#'
#' Note that this function will force the collection of \code{x}
#'
#' @export
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

