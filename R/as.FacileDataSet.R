# User-friendly functions used in the creation of a FacileDataSet
#
# The easiest way to create a FacileDataSet is to start from a well manicured
# SummarizeExperiment (or list of them).
#' Converts bioconductor assay containers into a FacileDataSet.
#'
#' @description
#' This function assumes you are only extracting one assay from the assay
#' container and creating a `FacileDataSet` from it. If the bioc-containers
#' you are using can support more than one assay, specify which one to extract
#' using the `source_assay` parameter, otherwise the first assay in the
#' container will be taken. The `source_assay` you are using to initialize the
#' `FacileDataSet` here will become its `default_assay` when its used later.
#'
#' The `feature_info` table is populated by the particular `fData` (`mcols`,
#' `y$genes`, etc.) from the assay container. The `pData` from the assay
#' container will be ingested as well, so do your best to ensure that this
#' represents the most complete version of the pData for the FacileDataSet
#' you will be creating.
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
#' @section Sample Covariates:
#'
#' The `pData` data.frame object will be picked off from all of the containers
#' provided in the list of datasets you are using to create the FacileDataSet.
#' `dataset` and `sample_id` columns will be forcibly added (or modified) as
#' columns to all of the individual `pData` data.frames.
#'
#' In order to insert the entirety of the `pData` elements into the inernal
#' `sample_covariate` table, we rely on the `dplyr::bind_rows` function to
#' create an uber `data.frame` which will be converted into an
#' entity-attribute-value table. Note that when row-binding, columns are matched
#' by name, and any missing columns with be filled with `NA`.
#'
#' Please ensure that the covariates across the `pData` data.frames have already
#' been harmonized!
#'
#' @section Feature meta-information:
#'
#' The feature information (aka "fData") are stored in an internal
#' `feature_info` SQLite table within the `FacileDataSet`. The information to
#' populate this table will be retrieved from the corresponding `fData`-like
#' `data.frame` from **the first** given bioc-container in the list.
#' This `data.frame` must define the following columns:
#'
#' * "feature_type": `character` (c('entrez', 'ensgid', 'enstid', 'genomic', 'custom'))
#' * "feature_id": `character`
#' * "name": `character`
#' * "meta": `character`
#' * "effective_length": `numeric`
#' * "source": `character`
#'
#' @md
#' @rdname as.FacileDataSet
#' @export
#'
#' @param x The bioconductor assay container to extract data from
#' @param path the directory to create the faciledataset into. Will create
#'   a default directory in the current working directory if not specified.
#'   This directory should not yet exist.
#' @param assay_name The name to use when storing the assay matrix from
#'   `x` into the faciledataset.
#' @param assay_type what type of assay is this? rnaseq, microarray, nanostring,
#'   isoseq (isoform expression), etc.
#' @param source_assay the name of the assay element in `x` to extract
#'   for use.
#' @param organism This is used to fetch the appropriate genesets when this
#'   dataset is used with the facileexplorer. Currently supported values are:
#'   `c("Homo sapiens", "Mus musculus", "unspecified")`.
#' @param dataset_name the `name` attribute of the FacileDataSet `meta.yaml`
#'   file.
#' @param ... more args
#' @return a [FacileDataSet()]
as.FacileDataSet <- function(x, path, assay_name, assay_type,
                             source_assay = NULL,
                             dataset_name = "DEFAULT_NAME",
                             organism = c("unspecified", "Homo sapiens", "Mus musculus"),
                             ...) {
  organism = match.arg(organism)
  UseMethod("as.FacileDataSet")
}

## These are the containers we can extract data from
## package=class
legit.as.classes <- c(
    "SummarizedExperiment"="SummarizedExperiment",
    "SummarizedExperiment"="RangedSummarizedExperiment",
    "Biobase"="ExpressionSet",
    "edgeR"="DGEList")

#' @method as.FacileDataSet default
#' @export
as.FacileDataSet.default <- function(x, ...) {
    xclass = class(x)[1L]
    if (! xclass %in% legit.as.classes) {
        stop("as.FacileDataSet not defined for object of class: ", xclass)
    }
    as.FacileDataSet(list(x))
}

#' @method as.FacileDataSet list
#' @export
#' @rdname as.FacileDataSet
as.FacileDataSet.list <- function(x, path, assay_name, assay_type,
                                  source_assay,
                                  organism,
                                  dataset_name, ...) {
  stopifnot(is.list(x))
  stopifnot(length(x) >= 1L)
  if (file.exists(path)) {
    stop("Target directory already exists: ", path)
  }

  # names(x) defines the name of the `dataset` value for each internal dataset
  if (length(x) > 1L) {
    # ensure that names of datasets are specified and unique
    stopifnot(
      is.character(names(x)),
      sum(duplicated(names(x))) == 0L)
  } else if (!is.character(names(x))) {
    # there's only one dataset here and no name is provided, so we use a
    # default dataset name if the list is unnamed
    names(x) <- tolower(fclass)
  }

  # All elements in list must be the same class, and a legit class at that!
  first <- x[[1L]]
  fclass <- class(first)[1L]
  # ... same class
  same.classes <- sapply(x, function(xx) class(xx)[1L] == fclass)
  stopifnot(all(same.classes))
  # ... legit class
  stopifnot(fclass %in% legit.as.classes)

  # load required namespace to deal with object of type `fclass`
  pkg <- local({
    pkg.idx <- which(sapply(legit.as.classes, function(clz) is(first, clz)))
    names(legit.as.classes)[pkg.idx]
  })
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " package required to convert these objects to a FacileDataSet")
  }

  # feature-space of assay containers must be identical
  same.fspace <- sapply(x, function(xx) {
    nrow(xx) == nrow(first) && all.equal(rownames(xx), rownames(first))
  })
  stopifnot(all(same.fspace))

  finfo <- fdata(first, validate = TRUE)
  pdats <- lapply(names(x), function(dname) {
      obj <- x[[dname]]
      out <- pdata(obj)
      out$dataset <- dname # Avoid dplyr::mutate to keep label attribute
      out$sample_id <- colnames(obj)
      dplyr::select(out, dataset, sample_id, everything())
  })

  pdat <- dplyr::bind_rows(pdats)
  col_descs = unlist(lapply(pdats, attr, which = "label"))
  col_descs = col_descs[!duplicated(names(col_descs))]
  attr(pdat, "label") <- col_descs

  eav.meta <- eav_metadata_create(pdat, covariate_def = NULL)

  adat <- lapply(names(x),
                 function(dname) {
                     adata(x[[dname]], assay = source_assay)
                 })

  ds_list = sapply(names(x), function(dname) {
      ds_annot(x[[dname]])
  }, simplify = FALSE)

  meta <- list(
    name = dataset_name,
    organism = organism,
    default_assay = assay_name,
    datasets = ds_list,
    sample_covariates = eav.meta)
  list(fdata = finfo, pdata = pdat, meta = meta, adata = adat)
}

#' Bioc-container specific data set annotation extraction functions
#'
#' Takes dataset-level annotion as stored by each type. DGEList has
#' no such slot, unfortunately, and thus gets the default. SE has a
#' metadata slot and can provide url and description. eSet just has
#' a character annotation and can provide a description.
ds_annot <- function(x, validate = FALSE, ...) {
  UseMethod("ds_annot")
}

ds_annot.SummarizedExperiment <- function(x, validate = FALSE) {
    out = metadata(x)
    if (validate && !(is.list(x) && setequal(names(x), c("url", "description")))) {
        stop("SummarizedExperiment data set level annotation was not a list with 'url' and 'description'.")
    }
    out
}

ds_annot.ExpressionSet <- function(x, validate = FALSE) {
    out = list(
        url = "http://www.google.com",
        description = annotation(x)
        )
    out
}

ds_annot.default <- function(x, validate = FALSE) {
    list(url = "http://google.com", description = "NoDescription")
}

#' Bioc-container specific fData extraction functions
#'
#' not exported on purpose
fdata <- function(x, validate = FALSE, ...) {
  UseMethod("fdata")
}
fdata.SummarizedExperiment <- function(x, validate = FALSE, ...) {
  stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE))
  out <- as.data.frame(SummarizedExperiment::mcols(x))
  if (validate) validate.fdata(out, ...) else out
}
fdata.ExpressionSet <- function(x, validate = FALSE, ...) {
  stopifnot(requireNamespace("Biobase", quietly = TRUE))
  out <- Biobase::fData(x)
  if (validate) validate.fdata(out, ...) else out
}
fdata.DGEList <- function(x, validate = FALSE, ...) {
  # stopifnot(requireNamespace("edgeR", quietly = TRUE))
  out <- x$genes
  if (validate) validate.fdata(out, ...) else out
}
validate.fdata <- function(x, ...) {
  stopifnot(is.data.frame(x))
  req.cols <- c(
    feature_type="character",
    feature_id="character",
    name="character",
    meta="character",
    effective_length="numeric")

  # Ensure required columns are there, and that they are of the right type
  for (cname in names(req.cols)) {
    vals <- x[[cname]]
    is.fn <- getFunction(paste0("is.", req.cols[cname]))
    if (is.null(vals)) {
      stop("Missing required fData column: ", cname)
    }
    if (!is.fn(vals)) {
      stop("'", cname, "' fData column not correct class: ", req.cols[cname])
    }
  }

  x[, names(req.cols)]
}

#' Bioc-container specific pData extraction functions
#'
#' not exported on purpose
pdata <- function(x, ...) {
  UseMethod("pdata")
}
pdata.SummarizedExperiment <- function(x, ...) {
    stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE))
    df = SummarizedExperiment::colData(x)
    ds = as.data.frame(df)
    if (length(metadata(df) > 0))
        attr(ds, "label") = unlist(metadata(df))
    validate.pdata(ds, ...)
}
pdata.ExpressionSet <- function(x, ...) {
  stopifnot(requireNamespace("Biobase", quietly = TRUE))
  validate.pdata(Biobase::pData(x), ...)
}
pdata.DGEList <- function(x, ...) {
  # stopifnot(requireNamespace("edgeR", quietly = TRUE))
  validate.pdata(x$samples, ...)
}
validate.pdata <- function(x, ...) {
  x
}

#' Bioc-container specific assay data extraction functions
#'
#' not exported on purpose
adata <- function(x, assay = NULL, ...) {
  UseMethod("adata")
}
adata.SummarizedExperiment <- function(x, assay = NULL, ...) {
  stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE))
  if (is.null(assay)) {
    assay <- SummarizedExperiment::assayNames(x)[1L]
  }
  if (!assay %in% SummarizedExperiment::assayNames(x)) {
    stop("Unknown assay in SummarizedExperiment: ", assay)
  }
  SummarizedExperiment::assay(x, assay)
}
adata.ExpressionSet <- function(x, assay = NULL, ...) {
  stopifnot(requireNamespace("Biobase", quietly = TRUE))
  if (is.null(assay)) {
    assay <- Biobase::assayDataElementNames(x)[1L]
  }
  if (!assay %in% Biobase::assayDataElementNames(x)) {
    stop("Unknown assay in ExpressoinSet: ", assay)
  }
  Biobase::assayDataElement(x, assay)
}
adata.DGEList <- function(x, assay = NULL, ...) {
  # stopifnot(requireNamespace("edgeR", quietly = TRUE))
  # DGEList only has one assay
  x$counts
}

# Internal functions to finalize as.FacileDataSet.* ============================

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
#'   complete description of what is expected in the `meta.yaml` file.
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
