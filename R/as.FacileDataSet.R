#' Converts bioconductor assay containers into a FacileDataSet.
#'
#' This function assumes you are only extracting one assay from the assay
#' container and creating a FacileDataSet from it. This requires that you
#' specify which assay (if the container has more than one) to extract, as
#' well as hand-crafting the feature_info correctly for the assay.
#'
#' @md
#' @rdname as.FacileDataSet
#' @export
#'
#' @param x The bioconductor assay container to extract data from
#' @param assay_name The name to use when storing the assay matrix from
#'   `x` into the faciledataset.
#' @param assay_type what type of assay is this? rnaseq, microarry, nanostring,
#'   isoseq (isoform expression), etc.
#' @param feature_info a `data.frame` that describes the information for the
#'   features (rows) of the assay you are extracting. Currently you had to
#'   hand-craft this. In the future we will provide automated default
#'   fData, rowData, etc. extractors for the source assay cointaner.
#' @param feature_type `c('entrez', 'ensgid', 'enstid', 'genomic', 'custom')`
#' @param metayaml a yaml file (or list of lists) that describes the covariates
#'   in the pData/colData of `x`. If not provided, a default one will be
#'   generated
#' @param organism This is used to fetch the appropriate genesets when this
#'   dataset is used with the facileexplorer. Currently supported values are:
#'   `c("Homo sapiens", "Mus musculus", "unspecified")`.
#' @param path the directory to create the faciledataset into. Will create
#'   a default directory in the current working directory if not specified.
#' @param source_assay the name of the assay element in `x` to extract
#'   for use.
#' @param ... more args
#' @return a [FacileDataSet()]
as.FacileDataSet <- function(x, assay_name, assay_type, feature_info,
                             feature_type, metayaml=NULL,
                             organism="unspecified",
                             path=tempfile("FacileDataSet-", getwd()),
                             source_assay=NULL, ...) {
  UseMethod("as.FacileDataSet")
}

##' @method as.FacileDataSet default
##' @export
as.FacileDataSet.default <- function(x, assay_name, assay_type, feature_info,
                                     feature_type, metayaml=NULL,
                                     organism="unspecified",
                                     path=tempfile("FacileDataSet-", getwd()),
                                     source_assay=NULL, ...) {
  stop("as.FacileDataSet not defined for object of class: ", class(x)[1L])
}

##' @method as.FacileDataSet ExpressionSet
##' @export
##' @rdname as.FacileDataSet
as.FacileDataSet.ExpressionSet <- function(x, assay_name, assay_type, feature_info,
                                           feature_type, metayaml=NULL,
                                           organism="unspecified",
                                           path=tempfile("FacileDataSet-", getwd()),
                                           source_assay=assayDataElementNames(x)[1L],
                                           ...) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package required to convert ExpresionSet to FacileDataSet")
  }
}

##' @method as.FacileDataSet SummarizedExperiment
##' @export
##' @rdname as.FacileDataSet
as.FacileDataSet.SummarizedExperiment <- function(x, assay_name, assay_type, feature_info,
                                                  feature_type, metayaml=NULL,
                                                  organism="unspecified",
                                                  path=tempfile("FacileDataSet-", getwd()),
                                                  source_assay=NULL, ...) {
}

##' @method as.FacileDataSet DGEList
##' @export
##' @rdname as.FacileDataSet
as.FacileDataSet.DGEList <- function(x, assay_name, assay_type, feature_info,
                                     feature_type, metayaml=NULL,
                                     organism="unspecified",
                                     path=tempfile("FacileDataSet-", getwd()),
                                     source_assay=NULL, ...) {
}

as.FacileDataSet.matrix <- function(x, assay_name, assay_type, feature_info,
                                    feature_type, metayaml=NULL,
                                    organism="unspecified",
                                    path=tempfile("FacileDataSet-", getwd()),
                                    source_assay=NULL, ...) {

}

##' @method as.FacileDataSet list
##' @export
##' @rdname as.FacileDataSet
as.FacileDataSet.list <- function(x, path, organism, assays=NULL, metayaml=NULL,
                                  ...) {
}

#' Create a facile covariate definition file from a sample `pData` data.frame
#'
#' Sample covariates (aka `pData`) are encoded in an
#' [entity-attribute-value (EAV) table](https://en.wikipedia.org/wiki/Entity-attribute-value_model).
#' Metadata about these covariates are stored in a `meta.yaml` file in the
#' `FacileDataSet` directory which enables the `FacileDataSet` to cast the value
#' stored in the EAV table to its native R type. This function generates the
#' list-of-list structure to represent the `sample_covariates` section of the
#' `meta.yaml` file.
#'
#' For simple `pData` covariates, each column is treated independantly from the
#' rest. There are some types of covariates which require multiple columns for
#' proper encoding, such as encoding survival information (see the **Survival**
#' section below). In these cases, the caller needs to provide an entry in the
#' `covariate_def` list that describes which `pData` columns (`varname`) goes
#' into the single facile covariate value. Please reference the
#' **Defining EAV Encodings** section for more information.
#'
#' keeping track of the "time to event" in one column, and a separae column to
#' indicate whether or not the event was a "death" or a censored. Still, these
#' data are stored in the single "value" column of the FacileDataSet's internal
#' entity-attribute-value (`sample_covariate`) table. In order to
#' encode these types of columns correctly, we need to provide more information
#' via the `covariate_def` parameter of this function.
#'
#' @section Survival:
#' Survival data is encoded by two columns. One column to indicate the
#' "time to event" and a second to indicate whether or not the denoted
#' tte is an "event" (1) or "censored" (0). The pair of columns will be encoded
#' into the `FacileDataSet`'s EAV table as a single (numeric) value. The
#' absolute value of the numeric indicates the "time to event" and the sign of
#' the value indicates its censoring status (If there are such data in `x`, it
#' must be in a (`tte_OS`, `event_OS`) pair of columns for "ordinary survival"
#' or a (`tte_PFS`, `event_PFS`) for progression free survival. Please see the
#' "Defining EAV Encodings" sections for more details.
#'
#' @md
#' @rdname eav-metadata
#' @export
#'
#' @param x a `pData` `data.frame`
#' @param covariate_def a named list of covariate definitions. The names of
#'   this list are the names the covariates will be called in the target
#'   `FacileDataSet`. The values of the list are:
#'   * `varname`: a `character()` of the column name(s) in `x` that this
#'      sample covariate was derived from. If more than one column is to be used
#'      for the facile covariate conversion (eg. if we are encoding survival),
#'      then provide a `length() > 1` character vector with the names of the
#'      columns in `x` that were used for the encoding. If this were encoding
#'      survival this might be `c("time", "event") columns, in that order.
#'   * label: a human readable label to use for this covariate in user facing
#'      scenarios in the facileverse.
#'   * `class`: the "facile class" of the covariate. This can either be
#'     `categorical`, `real`, or `right_censored` (for survival).
#'   * `levels`: (optional) if you want a `categorical` to be treated as a
#'     factor if it isn't already encoded as such in the `pData` itself, or if
#'     you want to rearrange the factor levels.
#'   * `type`: (optinal) this is used a a "grouping" level, particularly in
#'     the FacileExplorer. Not including this won't matter.
#'     TODO: talk about covariate groupings in
#'     `FacileExplorer::categoricalCovariateSelect`
#' @return a list-of-lists that encodes the `sample_covariate` section of the
#'   `meta.yaml` file for a `FacileDataSet`.
#'
#' @examples
#' # covariate_def definition to take tte_OS and tte_event columns and turn
#' # into a facile "OS" right_censored survival covariate
#' cc <- list(
#'   OS=list(
#'     varname=c("tte_OS", "tte_event"),
#'     label="Overall Survival",
#'     class="right_censored",
#'     type="clinical",
#'     description="Overall survival in days"))
create_eav_metadata <- function(x, covariate_def = list()) {
  stopifnot(is.data.frame(x))
  if (is.null(covariate_def)) covariate_def <- list()
  stopifnot(is.list(covariate_def))
  if (length(covariate_def)) {
    validate_covariate_def_list(covariate_def, x)
  }

  if ("dataset" %in% colnames(x)) x[['dataset']] <- NULL
  if ("sample_id" %in% colnames(x)) x[['sample_id']] <- NULL

  # generate generic covariate definitions for all columns
  gcd <- lapply(colnames(x), function(name) eavdef_for_column(x, name))
  names(gcd) <- colnames(x)

  # remove definitions in gcd that are provided in covariate_def
  axe <- lapply(covariate_def, '[[', 'varname')
  axe <- unique(unlist(axe, recursive = TRUE, use.names = FALSE))
  gcd[axe] <- NULL

  out <- c(gcd, covariate_def)
  out
}

#' @rdname eav-metadata
#'
#' @section Defining EAV Encodings:
#'
validate_covariate_def_list <- function(x, pdata) {
  # this is named a list of lists
  stopifnot(is.list(x), is.character(names(x)))
  is.lists <- sapply(x, is.list)
  stopifnot(all(is.lists))
  # the names are unique
  stopifnot(length(unique(names(x))) == length(x))

  # each list item has the following elements
  required.elements <- c('varname', 'class')
  # 'label', 'type' (group), and 'description' are also nice to have, but not
  # strictly required.
  has.elems <- sapply(x, function(el) all(required.elements %in% names(el)))
  stopifnot(all(has.elems))

  # varname entries are valid columns in `pdata`
  for (element in names(x)) {
    cnames <- x[[element]]$varname
    if (!is.character(cnames)) {
      stop(sprintf("%s:varname is not a character", element))
    }
    bad.cnames <- setdiff(cnames, colnames(pdata))
    if (length(bad.cnames)) {
      msg <- "%s:varname element '%s' not found as column in `pdata`"
      msg <- sprintf(msg, element, paste(bad.cnames, collapse="&"))
      stop(msg)
    }
  }
  TRUE
}

#' Generate entity-attribute-value definition for a column in a data.frame
#'
#' Creates the miniimal list-definition for a single column in a `pData`
#' `data.frame`. This function is not exported on purpose.
#'
#' @md
#' @param x a `data.frame`
#' @param column the name of the column in the `x` to parse.
#' @return a generic list-of-list defintion for `x[[column]]`
eavdef_for_column <- function(x, column) {
  vals <- x[[column]]
  if (is.null(vals)) stop("Unknown column in x: ", column)

  out <- list(varname=column, label=column, class='categorical', type="generic",
              description="no description provided")
  if (is.numeric(vals)) {
    out[['class']] <- "real"
  }
  if (is.factor(vals)) {
    out[['levels']] <- levels(vals)
  }
  out
}
