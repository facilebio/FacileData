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
#' create a `FacileDataSet` that spans multiple Bioc containters, i.e. you may
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
#' In order to insert the entirety of the `pData` elements into the internal
#' `sample_covariate` table, we rely on the `dplyr::bind_rows` function to
#' create an uber `data.frame` which will be converted into an
#' entity-attribute-value table. Note that when row-binding, columns are matched
#' by name, and any missing columns with be filled with `NA`.
#'
#' `ExpressionSet` pData `data.frames` should have an attribute called 'label', which
#' will be a named character vector with a description for each column. In the case of
#' a `SummarizedExperiment`, the `colData` should have named list in the `metadata`
#' slot with a character description of each column.
#'
#' `ExpressionSet`s should have a short textual description of the facet/dataset in
#' the `annotation` slot. Similarly, `SummarizedExperiment`s should have a list
#' in the `metadata` slot with `url` and `description` for the facet/dataset.
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
#' @param page_size parameter to tweak SQLite
#' @param cache_size parameter to tweak SQLite
#' @param chunk_rows parameter to tweak HDF5
#' @param chunk_cols parameter to tweak HDF5
#' @param chunk_compression parameter to tweak HDF5
#' @param ... more args
#' @return a [FacileDataSet()]
#' @importFrom yaml write_yaml
as.FacileDataSet <- function(x, path, assay_name, assay_type, source_assay,
                             dataset_name = "DEFAULT_NAME",
                             organism = c("unspecified", "Homo sapiens", "Mus musculus"),
                             page_size=2**12, cache_size=2e5,
                             chunk_rows=5000, chunk_cols="ncol",
                             chunk_compression=5,
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
                                  dataset_name,
                                  organism = c("unspecified", "Homo sapiens", "Mus musculus"),
                                  page_size=2**12, cache_size=2e5,
                                  chunk_rows=5000, chunk_cols="ncol",
                                  chunk_compression=5,
                                  ...) {
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
  same.fspace <- vapply(x, function(xx) {
    nrow(xx) == nrow(first) && all.equal(rownames(xx), rownames(first))
  }, logical(1))
  stopifnot(all(same.fspace))

  finfo <- fdata(first, validate = TRUE)

  pdats <- lapply(names(x), function(dname) {
      obj <- x[[dname]]
      out <- pdata(obj)
      out$dataset <- dname
      out$sample_id <- colnames(obj)
      dplyr::select(out, dataset, sample_id, everything())
  })

  pdat <- bind_pdata_rows(pdats)
  pdat$sample_id = gsub("__", "",  pdat$sample_id) # __ has as special meaning for Facile

  col_descs_list = lapply(x, pdata_metadata)
  col_descs = unlist(col_descs_list, FALSE, FALSE)
  names(col_descs) = unlist(sapply(col_descs_list, names, simplify = FALSE), FALSE, FALSE)
  col_descs = col_descs[!duplicated(names(col_descs))]

  eav.meta <- eav_metadata_merge(pdat, col_descs)

  adat <- sapply(names(x),
                 function(dname) {
                     ds = adata(x[[dname]], assay = source_assay)
                     colnames(ds) = gsub("__", "", colnames(ds)) # __ has as special meaning for Facile
                     ds
                 }, simplify = FALSE)


  ## Make YAML and Initialize FDS
  ds_list = sapply(names(x), function(dname) {
      ds_annot(x[[dname]])
  }, simplify = FALSE)

  meta <- list(
    name = dataset_name,
    organism = organism,
    default_assay = assay_name,
    datasets = ds_list,
    sample_covariates = col_descs
   )

  meta_yaml = paste0(tempfile(), ".yaml")
  write_yaml(meta, meta_yaml)
  fds <- initializeFacileDataSet(path, meta_yaml)

  ## FIXME: make pdat_eav
  pdat_eav = as.EAVtable(pdat, eav.meta)

  ## Check for duplicates as SQLite will raise exception
  stopifnot(nrow(duplicated(pdat_eav %>% select(dataset, sample_id, variable))) == 0)

  ## Register sample covariate info into Facile SQLite
  sample.covs <- pdat_eav %>%
      mutate(
          type = "general",  ## FIXME: care about type later
          date_entered = as.integer(Sys.time())
      ) %>% append_facile_table(fds, 'sample_covariate')

  ## Register sample info into Facile SQLite
  sample.info <- pdat_eav %>%
      select(dataset, sample_id) %>%
      distinct(dataset, sample_id) %>%
      mutate(parent_id =  "") %>%
        append_facile_table(fds, 'sample_info')

  ## insert the first assay
  if (is.integer(adat[[1]]))
      storage_mode = "integer"
  else
      storage_mode = "numeric"

  tstart <- Sys.time()
  samples <- addFacileAssaySet(
    fds,
    adat,
    facile_assay_name=meta$default_assay,
    facile_assay_type='rnaseq',
    facile_feature_type=finfo$feature_type[1],
    facile_feature_info=finfo,
    facile_assay_description=meta$name,
    storage_mode=storage_mode,
    chunk_rows=chunk_rows,
    chunk_cols=chunk_cols,
    chunk_compression=chunk_compression
  )
  tend <- Sys.time()
  message("Time taken: ", tend - tstart) ## ~ 30 seconds

  fds
}

#' Bioc-container specific data set annotation extraction functions
#'
#' Takes dataset-level annotion as stored by each type. DGEList has
#' no such slot, unfortunately, and thus gets the default. SE has a
#' metadata slot and can provide url and description. eSet just has
#' a character annotation and can provide a description.
#' @param x SummarizedExperiment, ExpressionSet or DGEList
#' @param validate single logical, check results
#' @param ... additional args (ignored for now)
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

#' BioC-container specific fData extraction functions
#'
#' not for export
#' @param x SummarizedExperiment, ExpressionSet or DGEList
#' @param validate single logical, check results
#' @param ... additional args (ignored for now)
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
    effective_length="numeric",
    source="character"
  )

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
#' not for export
#' @param x SummarizedExperiment, ExpressionSet or DGEList
#' @param ... additional args, ignored for now
pdata <- function(x, ...) {
  UseMethod("pdata")
}
pdata.SummarizedExperiment <- function(x, ...) {
    stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE))
    df = SummarizedExperiment::colData(x)
    ds = as.data.frame(df)
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

#' Bioc-container specific pData extraction functions
#'
#' Get metadata on columns of sample info data.frame (label, etc.) for
#' inclusion in metadata YAML.
#' not for export
#' @param x SummarizedExperiment, ExpressionSet or DGEList
#' @param ... additional args, ignored for now
pdata_metadata <- function(x, ...) {
  UseMethod("pdata_metadata")
}
pdata_metadata.SummarizedExperiment <- function(x, ...) {
    stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE))
    sinfo = SummarizedExperiment::colData(x)
    defs = S4Vectors::metadata(sinfo)
    defs
}
pdata_metadata.ExpressionSet <- function(x, ...) {
    sinfo = pdata(x)
    defs = attributes(sinfo)$label
    names(defs) = colnames(sinfo)
    defs = lapply(defs, function(el) { list(description = el) })
    defs
}
pdata_metadata.DGEList <- function(x, ...) {
    sinfo = x$samples
    defs = sapply(colnames(sinfo), function(el) {
        list(description = el, label = el, type = "general")
    }, simplify = FALSE)
    defs
}

#' Bioc-container specific assay data extraction functions
#'
#' Get assay matrix
#' @param x SummarizedExperiment, ExpressionSet or DGEList
#' @param ... additional args, ignored for now
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
