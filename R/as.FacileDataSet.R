#' Converts bioconductor assay containers into a FacileDataSet.
#'
#' @description
#' This function assumes you are only extracting one assay from the assay
#' container and creating a `FacileDataSet` from it. This requires that you
#' specify which assay (if the container has more than one) to extract. The
#' `feature_info` table is populated by the particular `fData` (`mcols`,
#' `y$genes`, etc.) from the assay container. The `pData` from the assay
#' container will be ingestged as well.
#'
#' **Please read the Details section**  for more a more complete description
#' of the data that is required to create a complete `FacileDataSet`
#'
#' @details
#'
#' A `FacileDataSet` can be created from a number of different Bioconductor
#' containers, such as a [Biobase::ExpressionSet],
#' [SummarizedExperiment::SummarizedExperiment], or an [edgeR::DGEList]. To
#' create a `FacileDataSet` that spans multiple Bioc containters, ie. you may
#' have one ExpressionSet per indication in the TCGA. You can make
#' `FacileDataSet` to encompass the data from *all* of these indications by
#' providing a `list` of `ExpressionSet`s. The `list` should have its `names()`
#' set to each of the TCGA indications ("BLCA", "BRCA", etc.) the data came
#' from.
#'
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

#' @method as.FacileDataSet default
#' @export
as.FacileDataSet.default <- function(x, assay_name, assay_type, feature_info,
                                     feature_type, metayaml=NULL,
                                     organism="unspecified",
                                     path=tempfile("FacileDataSet-", getwd()),
                                     source_assay=NULL, ...) {
  stop("as.FacileDataSet not defined for object of class: ", class(x)[1L])
}

#' @method as.FacileDataSet ExpressionSet
#' @export
#' @rdname as.FacileDataSet
as.FacileDataSet.ExpressionSet <- function(x, assay_name, assay_type, feature_info,
                                           feature_type, metayaml=NULL,
                                           organism="unspecified",
                                           path=tempfile("FacileDataSet-", getwd()),
                                           source_assay=assayDataElementNames(x)[1L],
                                           dset="eset", ...) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package required to convert ExpresionSet to FacileDataSet")
  }
  as.FacileDataSet(list(x), ...)
}

#' @method as.FacileDataSet SummarizedExperiment
#' @export
#' @rdname as.FacileDataSet
as.FacileDataSet.SummarizedExperiment <- function(x, assay_name, assay_type, feature_info,
                                                  feature_type, metayaml=NULL,
                                                  organism="unspecified",
                                                  path=tempfile("FacileDataSet-", getwd()),
                                                  source_assay=NULL,
                                                  dset="sumexp", ...) {
  as.FacileDataSet(list(x), ...)
}

#' @method as.FacileDataSet DGEList
#' @export
#' @rdname as.FacileDataSet
as.FacileDataSet.DGEList <- function(x, assay_name, assay_type, feature_info,
                                     feature_type, metayaml=NULL,
                                     organism="unspecified",
                                     path=tempfile("FacileDataSet-", getwd()),
                                     source_assay=NULL, dset="dgelist", ...) {
  as.FacileDataSet(list(x), ...)
}

# These are the containers we can extract data from
legit.as.classes <- c(
  "SummarizedExperiment"="SummarizedExperiment",
  "Biobase"="ExpressionSet",
  "edgeR"="DGEList")

#' @method as.FacileDataSet list
#' @export
#' @rdname as.FacileDataSet
as.FacileDataSet.list <- function(x, path, organism, assays=NULL, metayaml=NULL,
                                  ...) {
  stopifnot(length(x) >= 1L)
  first <- x[[1L]]
  fclass <- class(first)[1L]

  # All elements in list must be the same class, and a legit class at that!
  # ... same class
  same.classes <- sapply(x, function(xx) class(xx)[1L] == fclass)
  stopifnot(all(same.classes))
  # ... legit class
  stopifnot(fclass %in% legit.as.classes)

  # load required namespace to deal with object of type `fclass`
  pkg.idx <- which(sapply(legit.as.classes, function(clz) is(first, clz)))
  pkg <- names(legit.as.classes)[pkg.idx]
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " package required to convert ExpresionSet to FacileDataSet")
  }

  # names(x) defines the name of the `dataset` value for each, uhm ... dataset
  if (length(x) > 1L) {
    # ensure that names of datasets are specified and unique
    stopifnot(
      is.character(names(x)),
      sum(duplicated(names(x))) == 0L)
  } else if (!is.character(names(x))) {
    # use a default dataset name if the list is unnamed
    names(x) <- tolower(fclass)
  }

  # feature-space of assay containers must be identical
  same.fspace <- sapply(x, function(xx) {
    nrow(xx) == nrow(first) && all.equal(rownames(xx), rownames(first))
  })
  stopifnot(all(same.fspace))

  # extract feature_info data.frame
  if (is(first, "SummarizedExperiment")) {
    finfo <- as.data.frame(SummarizedExperiment::mcols(first))
  } else if (is(first, "ExpressionSet")) {
    finfo <- Biobase::fData(first)
  } else if (is(first, "DGEList")) {
    finfo <- first$genes
  } else {
    stop("Shouldn't be here")
  }

}

# Internal functions to finalzie as.FacileDataSet.* ============================

#' Finalizes conversion of bioconductor containers into a FacileDataSet.
#'
#' This is the final call that the various
#' `as.FacileDataSet.(SummarizedExperiment)` functions delegate to in order to
#' actually create the `FacileDataSet` on disk. It is used to put all of the
#' minimal "bits" for a valid `FacileDataSet` into place, and is only meant to
#' initialize it with data from one assay type.
#'
#' To add more assays to the `FacileDataSet`, use the [addFacileAssaySet()]
#' function.
#'
#' **This function is intentionally not exported**, however a savvy user may
#' find themselves calling it to fill out a complete FacileDataSet after it has
#' been initially constructed via a call to `as.FacileDataSet(stuff, ...)`
#' not exported
#'
#' @md
#' @export
#'
#' @param path the path to the folder that will contain the facile contents
#' @param meta_file the path on disk where the `meta.yaml` file exists for
#'   this `FacileDataSet`. Reference the help in [FacileDataSet()] for a more
#'   complete description fo what is expected in the `meta.yaml` file.
#' @param page_size parameter to tweak SQLite
#' @param cache_size parameter to tweak SQLite
#' @param sample_covariates (named) list of pData `data.frame`s that will be
#'   inserted into the `FacileDataSet`.
#' @param assays a (named) list of assay matrices for the samples across the
#'   datasets
#' @param assay_name the name to use for the assay data in `assays`.
#' @param assay_type the type of assay (ie. `"rnaseq"`, `"affymetrix"`, etc.)
#' @param assay_feature_type the name of the "feature space" over the rows
#'   of this assay, something like `"entrez"`, `"ensgid"`, `"enstid"`, etc.
#' @param assay_feature_info the `feature_info` `data.frame` for the features
#'   in the rows of `assays`
as_FacileDataSet <- function(sample_covariates, assays, assay_name, assay_type,
                             assay_feature_type, assay_feature_info,
                             assay_description=assay_name,
                             path, meta_file, page_size=2**12, cache_size=2e5,
                             chunk_rows=5000, chunk_cols="ncol",
                             chunk_compression=5, covariate_def=list(),
                             metayaml=NULL, ...) {
  # combine list of `pData` data.frames into one, generate default yaml file
  # and generate long-form sample_covariate_table
  scovs <- df2eav(sample_covariates, covariate_def, metayaml, ...)

  # initialize directory structure with metat.yaml
  fds <- initializeFacileDataSet(path, meta_file, page_size, cache_size)

  # insert the first assay
  tstart <- Sys.time()
  samples <- addFacileAssaySet(
    fds,
    assays,
    facile_assay_name='rnaseq',
    facile_assay_type='rnaseq',
    facile_feature_type='entrez',
    facile_feature_info=facile_feature_info,
    facile_assay_description="Recounted RNA-seq data from first batch of Atezo Trials",
    storage_mode='integer',
    chunk_rows=5000, chunk_cols='ncol', chunk_compression=4)
  tend <- Sys.time()
  message("Time taken: ", tend - tstart) ## ~ 30 seconds
}

