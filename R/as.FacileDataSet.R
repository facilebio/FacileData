#' Converts a (list of) bioconductor assay containers into a FacileDataSet.
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
#' containers, such as a `Biobase::ExpressionSet`,
#' `SummarizedExperiment::SummarizedExperiment`, or an `edgeR::DGEList`. To
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
#' * "feature_type": `string`, one of: `"entrez"`, `"ensgid"`, `"enstid"`,
#'     `"genomic"`, `"custom"`.
#' * "feature_id": `string`
#' * "name": `string`
#' * "meta": `string`
#' * "effective_length": `integer`
#' * "source": `string`
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
#'   dataset is used with the facileexplorer. Put species name here, ie.
#'   `"Homo sapiens"`, `"Mus musculus"`, etc. Default is `"unspecified"`,
#'   which isn't really helpful.
#' @param dataset_name the `name` attribute of the FacileDataSet `meta.yaml`
#'   file.
#' @param dataset_meta a named (by names(x)) with meta data about the datasets
#'   that appear in the list of datasets `x`. List elements per dataset should
#'   minimally include a `description` and `url` string.
#' @param page_size parameter to tweak SQLite
#' @param cache_size parameter to tweak SQLite
#' @param chunk_rows parameter to tweak HDF5
#' @param chunk_cols parameter to tweak HDF5
#' @param chunk_compression parameter to tweak HDF5
#' @param ... more args
#' @return a [FacileDataSet()]
#' @importFrom yaml write_yaml
as.FacileDataSet <- function(x, path, assay_name, assay_type, source_assay,
                             assay_description = paste("Description for ", assay_name),
                             dataset_name = "DEFAULT_NAME",
                             dataset_meta = NULL,
                             organism = "unspecified",
                             page_size=2**12, cache_size=2e5,
                             chunk_rows=5000, chunk_cols="ncol",
                             chunk_compression=5, covariate_def = NULL,...) {
  UseMethod("as.FacileDataSet")
}

## These are the containers we can extract data from
## package=class

legit.as.classes <- tribble(
  ~package,               ~class,
  "SummarizedExperiment", "SummarizedExperiment",
  "SummarizedExperiment", "RangedSummarizedExperiment",
  "Biobase",              "ExpressionSet",
  "edgeR",                "DGEList",
  "limma",                "EList")

#' @method as.FacileDataSet default
#' @export
as.FacileDataSet.default <- function(x, ...) {
    xclass = class(x)[1L]
    if (!xclass %in% legit.as.classes[["class"]]) {
        stop("as.FacileDataSet not defined for object of class: ", xclass)
    }
    as.FacileDataSet(list(x), ...)
}

#' @method as.FacileDataSet list
#' @export
#' @rdname as.FacileDataSet
as.FacileDataSet.list <- function(x, path, assay_name, assay_type,
                                  source_assay = NULL,
                                  assay_description = paste("Description for ", assay_name),
                                  dataset_name = "DEFAULT_NAME",
                                  dataset_meta = list(),
                                  organism = "unspecified",
                                  page_size=2**12, cache_size=2e5,
                                  chunk_rows=5000, chunk_cols="ncol",
                                  chunk_compression=5,
                                  covariate_def = NULL, ...) {
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
    names(x) <- "dataset"
  }


  # All elements in list must be the same class, and a legit class at that!
  first <- x[[1L]]
  fclass <- class(first)[1L]
  # ... same class
  same.classes <- sapply(x, function(xx) class(xx)[1L] == fclass)
  stopifnot(all(same.classes))
  # ... legit class
  stopifnot(fclass %in% legit.as.classes[["class"]])

  # load required namespace to deal with object of type `fclass`
  pkg <- local({
    info <- filter(legit.as.classes, class == fclass)
    info[["package"]]
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

  # It's possible that expression containers have colnames(x) == NULL (ugh),
  # so let's provide default sample names if they are missing
  x <- lapply(x, function(xdat) {
    cnames <- colnames(xdat)
    if (is.null(cnames) || any(duplicated(cnames))) {
      colnames(xdat) <- paste0("sample_", seq(ncol(xdat)))
    }
    xdat
  })

  pdats <- lapply(names(x), function(dname) {
    obj <- x[[dname]]
    pdat <- pdata(obj)
    for (cname in colnames(pdat)) {
      vals <- pdat[[cname]]
      if (is(vals, "Surv")) pdat[[cname]] <- as_cSurv(vals)
    }
    pdat %>%
      mutate(dataset = dname, sample_id = colnames(obj)) %>%
      select(dataset, sample_id, everything())
  })

  pdat <- bind_rows(pdats)
  # __ has as special meaning for Facile
  pdat$sample_id = gsub("__", "",  pdat$sample_id)

  # col_descs_list = lapply(x, pdata_metadata, covariate_metadata)
  # col_descs = unlist(col_descs_list, FALSE, FALSE)
  # names(col_descs) = unlist(sapply(col_descs_list, names, simplify = FALSE),
  #                           recursive = FALSE, use.names = FALSE)
  # col_descs = col_descs[!duplicated(names(col_descs))]
  #
  # eav.meta <- eav_metadata_merge(pdat, col_descs)

  eav = as.EAVtable(pdat, covariate_def = covariate_def)
  eav_meta <- attr(eav, "covariate_def")

  # Check for duplicates entries in the eav table here, otherwise injecting it
  # into the database will raise an error
  stopifnot(
    nrow(distinct(eav, dataset, sample_id, variable)) == nrow(eav))

  # list of assay matricess
  adat <- sapply(names(x), function(dname) {
    ds = adata(x[[dname]], assay = source_assay)
    # __ has as special meaning for Facile
    colnames(ds) = gsub("__", "", colnames(ds))
    ds
  }, simplify = FALSE)


  # Make YAML and Initialize FDS
  ds_list <- sapply(names(x), function(dname) {
    ds_annot(x[[dname]], dataset_meta[[dname]])
  }, simplify = FALSE)

  meta <- list(
    name = dataset_name,
    organism = organism,
    default_assay = assay_name,
    datasets = ds_list,
    sample_covariates = eav_meta)

  meta_yaml <- paste0(tempfile(), ".yaml")
  yaml::write_yaml(meta, meta_yaml)
  path <- initializeFacileDataSet(path, meta_yaml)

  assert_directory_exists(path, access = "w")
  fds <- FacileDataSet(path)

  # add sample covariates to table
  sample.covs <- eav %>%
    mutate(date_entered = as.integer(Sys.time())) %>%
    append_facile_table(fds, "sample_covariate")

  # Add samples to sample_info table
  sample.info <- eav %>%
    distinct(dataset, sample_id) %>%
    mutate(parent_id = "") %>%
    append_facile_table(fds, "sample_info")

  # insert the first assay
  if (is.integer(adat[[1]])) {
    storage_mode = "integer"
  } else {
    storage_mode = "numeric"
  }

  tstart <- Sys.time()
  samples <- addFacileAssaySet(
    fds,
    adat,
    facile_assay_name=meta$default_assay,
    facile_assay_type=assay_type,
    facile_feature_type=finfo$feature_type[1],
    facile_feature_info=finfo,
    facile_assay_description=assay_description,
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
#' @param meta a list of description stuff for the dataset, this can act to
#'   override what's there, already
#' @param ... additional args (ignored for now)
ds_annot <- function(x, meta = NULL, validate = FALSE, ...) {
  UseMethod("ds_annot")
}

ds_annot.SummarizedExperiment <- function(x, meta = NULL, validate = FALSE,
                                          ...) {
  ns4 <- tryCatch(loadNamespace("S4Vectors"), error = function(e) NULL)
  if (is.null(ns4)) stop("S4Vectors required")

  out <- .merge_ds_annot(ns4$metadata(x), meta)
  if (validate) .ds_annot_validate(out)
  out
}

ds_annot.ExpressionSet <- function(x, meta = NULL, validate = FALSE, ...) {
  out <- .merge_ds_annot(NULL, meta)
  if (validate) .ds_annot_validate(out)
  out
}

ds_annot.default <- function(x, meta = NULL, validate = FALSE, ...) {
  out <- .merge_ds_annot(NULL, meta)
  if (validate) .ds_annot_validate(out)
  out
}

.merge_ds_annot <- function(x, y) {
  if (is.null(x)) x <- list()
  if (is.null(y)) y <- list()
  elems <- union(names(x), names(y))
  out <- sapply(elems, function(name) {
    c(x[[name]], y[[name]])[1L]
  }, simplify = FALSE)
  if (is.null(out[["description"]])) {
    out[["description"]] <- "No Description Provided"
  }
  if (is.null(out[["url"]])) {
    out[["url"]] <- "http://example.com"
  }
  out
}

.ds_annot_validate <- function(x) {
  if (!is.list(x) && test_subset(c("url", "description"), names(x))) {
    stop("Dataset metadata needs to be a list, minimally with 'url' and ",
         "'description' components")
  }
  invisible(x)
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
  ns <- tryCatch(loadNamespace("SummarizedExperiment"), error = function(e) NULL)
  if (is.null(ns)) stop("SummarizedExperiment required")
  ns4 <- tryCatch(loadNamespace("S4Vectors"), error = function(e) NULL)
  if (is.null(ns4)) stop("S4Vectors required")

  # out <- ns$as.data.frame(ns$mcols(x))
  out <- ns4$as.data.frame.DataTable(ns$rowData(x))
  if (validate) validate.fdata(out, ...) else out
}
fdata.ExpressionSet <- function(x, validate = FALSE, ...) {
  ns <- tryCatch(loadNamespace("Biobase"), error = function(e) NULL)
  if (is.null(ns)) stop("Biobase required")

  out <- ns$fData(x)
  if (validate) validate.fdata(out, ...) else out
}
fdata.DGEList <- function(x, validate = FALSE, ...) {
  # stopifnot(requireNamespace("edgeR", quietly = TRUE))
  out <- x$genes
  if (validate) validate.fdata(out, ...) else out
}
fdata.EList <- function(x, validate = FALSE, ...) {
  out <- x$genes
  if (validate) validate.fdata(out, ...) else out
}

validate.fdata <- function(x, ...) {
  stopifnot(is.data.frame(x))
  if (is.null(x[["source"]])) {
    x[["source"]] <- rep("unspecified", nrow(x))
  }

  req.cols <- c(
    feature_type="character",
    feature_id="character",
    name="character",
    meta="character",
    # effective_length="numeric",
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
pdata <- function(x, covariate_metadata = NULL, ...) {
  UseMethod("pdata")
}

#' @noRd
pdata.SummarizedExperiment <- function(x, covariate_metadata = NULL,  ...) {
  ns <- tryCatch(loadNamespace("SummarizedExperiment"), error = function(e) NULL)
  if (is.null(ns)) stop("SummarizedExperiment required")
  ns4 <- tryCatch(loadNamespace("S4Vectors"), error = function(e) NULL)
  if (is.null(ns4)) stop("S4Vectors required")

  df <- ns$colData(x)
  ds <- ns4$as.data.frame.DataTable(df)
  validate.pdata(ds, ...)
}

#' @noRd
pdata.ExpressionSet <- function(x, covariate_metadata = NULL,  ...) {
  ns <- tryCatch(loadNamespace("Biobase"), error = function(e) NULL)
  if (is.null(ns)) stop("Biobase required")
  validate.pdata(ns$pData(x), ...)
}

#' @noRd
pdata.DGEList <- function(x, covariate_metadata = NULL,  ...) {
  # stopifnot(requireNamespace("edgeR", quietly = TRUE))
  ignore.cols <- c("lib.size", "norm.factors")
  grp <- x$samples$group
  if (is.null(grp) || any(is.na(grp)) || all(grp == 1)) {
    ignore.cols <- c(ignore.cols, "group")
  }
  validate.pdata(x$samples[, !colnames(x$samples) %in% ignore.cols], ...)
}

#' @noRd
pdata.EList <- function(x, covariate_metadata = NULL,  ...) {
  validate.pdata(x$targets, ...)
}

validate.pdata <- function(x, ...) {
  x
}

#' Bioc-container specific pData extraction functions
#'
#' Get metadata on columns of sample info data.frame (label, etc.) for
#' inclusion in metadata YAML.
#'
#' @param x SummarizedExperiment, ExpressionSet or DGEList
#' @param ... additional args, ignored for now
#' @export
pdata_metadata <- function(x, ...) {
  UseMethod("pdata_metadata")
}

#' SummarizedExperiment method
#' @export
#' @noRd
pdata_metadata.SummarizedExperiment <- function(x, ...) {
  ns <- tryCatch(loadNamespace("SummarizedExperiment"), error = function(e) NULL)
  if (is.null(ns)) stop("SummarizedExperiment required")
  ns4 <- tryCatch(loadNamespace("S4Vectors"), error = function(e) NULL)
  if (is.null(ns4)) stop("S4Vectors required")

  sinfo <- ns$colData(x)
  ns4$metadata(sinfo)
}

#' ExpressionSet method
#' @export
#' @noRd
pdata_metadata.ExpressionSet <- function(x, ...) {
  sinfo <- pdata(x)
  defs <- attributes(sinfo)$label
  names(defs) <- colnames(sinfo)
  defs <- lapply(defs, function(el) { list(description = el) })
  defs
}

#' @export
#' @noRd
pdata_metadata.DGEList <- function(x, ...) {
  sinfo <- x$samples
  defs <- sapply(colnames(sinfo), function(el) {
    list(description = el, label = el, type = "general")
  }, simplify = FALSE)
  defs
}

#' @export
#' @noRd
pdata_metadata.EList <- function(x, ...) {
  sinfo <- x$targets
  defs <- sapply(colnames(sinfo), function(el) {
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

#' @noRd
adata.SummarizedExperiment <- function(x, assay = NULL, ...) {
  ns <- tryCatch(loadNamespace("SummarizedExperiment"), error = function(e) NULL)
  if (is.null(ns)) stop("SummarizedExperiment required")
  anames <- ns$assayNames(x)
  if (is.null(assay)) assay <- anames[1L]
  if (!assay %in% anames) {
    stop("Unknown assay in SummarizedExperiment: ", assay)
  }
  ns$assay(x, assay)
}

#' @noRd
adata.ExpressionSet <- function(x, assay = NULL, ...) {
  ns <- tryCatch(loadNamespace("Biobase"), error = function(e) NULL)
  if (is.null(ns)) stop("Biobase required")
  anames <- ns$assayDataElementNames(x)
  if (is.null(assay)) assay <- anames[1L]
  if (!assay %in% anames) {
    stop("Unknown assay in ExpressoinSet: ", assay)
  }
  ns$assayDataElement(x, assay)
}

#' @noRd
adata.DGEList <- function(x, assay = NULL, ...) {
  # stopifnot(requireNamespace("edgeR", quietly = TRUE))
  # DGEList only has one assay
  x$counts
}

#' @noRd
adata.EList <- function(x, assay = NULL, ...) {
  # EList only has one assay
  x[["E"]]
}
