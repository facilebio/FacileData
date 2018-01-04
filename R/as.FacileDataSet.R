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
