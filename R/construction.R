## The code in here builds FacileDataSets from Bioconductor like objects

#' Create an empty FacileDataSet
#'
#' This is a helper function that is currently only called from
#' `as.FacileDataSet`
#'
#' @importFrom rhdf5 h5createFile h5createGroup H5close
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery
#' @importFrom RSQLite SQLite
#'
#'
#' @param path the directory to create which will house the
#'   \code{FacileDataSet}
#' @param covariate_definition the path to the covariate definition file
#' @param page_size,cache_size \code{pragma} values to setup the backend SQLite
#'   database
#' @return inivisibly returns the path to the successfully created datastore
initializeFacileDataSet <- function(path, meta_file,
                                    page_size=2**12, cache_size=2e5) {
  assert_valid_meta_file(meta_file)
  path <- normalizePath(path, mustWork=FALSE)
  if (file.exists(path) || dir.exists(path)) {
    stop("Collision with file: ", path)
  }
  parent.dir <- dirname(path)
  assert_directory_exists(parent.dir)
  assert_access(parent.dir, 'w') ## check parent directory is writeable

  dir.create(path)
  dir.create(file.path(path, 'custom-annotation'))
  # create dummy file in custom-annotation directory for source control and
  # s3 synching of empty directories
  writeLines("Empty placeholder for synching empty directories",
             file.path(path, 'custom-annotation', 'README.txt'))
  ## file.copy(covariate_definition, file.path(path, 'sample-covariate-info.yaml'))
  file.copy(meta_file, file.path(path, 'meta.yaml'))

  ## Create sqlite db
  db.fn <- file.path(path, 'data.sqlite')
  con <- dbConnect(SQLite(), db.fn)
  on.exit(dbDisconnect(con))
  dbExecute(con, 'pragma temp_store=MEMORY;')
  dbExecute(con, sprintf('pragma page_size=%d', page_size))
  dbExecute(con, sprintf('pragma cache_size=%d;', cache_size))
  sql.fn <- system.file('extdata', 'init', 'faciledataset.sql',
                        package='FacileData')
  db.sql <- sqlFromFile(sql.fn)
  executeSQL(con, db.sql)

  ## Create empty HDF5 file
  hd5.fn <- file.path(path, 'data.h5')
  h5createFile(hd5.fn)
  h5createGroup(hd5.fn, 'assay')
  invisible(path)
}

#' @export
#' @importFrom tools file_ext
assert_valid_meta_file <- function(fn, as.list = FALSE) {
  assert_file(fn)
  if (!tolower(file_ext(fn)) ==  'yaml') {
    stop("meta file must be a yaml file")
  }
  dat <- yaml.load_file(fn)
  req.toplevel <- c('name', 'organism', 'datasets', 'sample_covariates',
                    'default_assay')
  miss.toplevel <- setdiff(req.toplevel, names(dat))
  if (length(miss.toplevel)) {
    stop("Missing the following definitions in meta file: ",
         paste(miss.toplevel, collapse=","))
  }
  if (as.list) dat else fn
}

.feature.types <- c('entrez', 'ensgid', 'enstid', 'custom')
.assay.types <- c(
  "rnaseq", "pseudobulk", "isoseq", "normcounts", "nanostring",
  "qpcrct", "qpcrdct", "lognorm", "raw")

.storage.modes <- c('integer', 'numeric')

supported.assay.container <- function(x) {
  supported <- c("DGEList", "EList", "eSet", "SummarizedExperiment", "matrix")
  any(is(x) %in% supported)
}

extract.assay <- function(x, assay_name=NULL) {
  if (is(x, 'DGEList')) {
    out <- x$counts
  } else if (is(x, "EList")) {
    out <- x$E
  } else if (is(x, 'eSet')) {
    ns <- loadNamespace("Biobase")
    if (is.null(assay_name)) assay_name <- assayDataElementNames(x)[1]
    out <- ns$assayDataElement(x, assay_name)
  } else if (is(x, 'SummarizedExperiment')) {
    ns <- loadNamespace("SummarizedExperiment")
    out <- if (is.null(assay_name)) {
      ns$assays(x)[[1L]]
    } else {
      ns$assay(x, assay_name)
    }
  } else if (is(x, 'matrix')) {
    out <- x
  } else {
    stop("Unsupported assay container: ", class(x)[1L])
  }
  stopifnot(is(out, 'matrix'),
            nrow(out) == nrow(x),
            ncol(out) == ncol(x),
            all.equal(rownames(out), rownames(x)),
            all.equal(colnames(out), colnames(x)))
  out
}

assert_valid_assay_datasets <- function(datasets, facile_feature_info,
                                        storage_mode, assay_name=NULL) {
  supported <- sapply(datasets, supported.assay.container)
  if (any(!supported)) {
    idxs <- which(!supported)
    bad <- sapply(datasets[idxs], function(d) class(d)[1L]) |> unique()
    stop("Unsupported assay container(s):", paste(bad, collapse=","))
  }
  dnames <- names(datasets)
  if (is.null(dnames) ||
      all(dnames != make.names(dnames)) ||
      length(unique(dnames)) != length(dnames)) {
    # stop("names(datasets) must be valid varnames and all unique")
    warning("names(datasets) are not valid R variable names", immediate. = TRUE)
  }

  smodes <- sapply(datasets, function(d) {
    class(extract.assay(d, assay_name)[1L])
  })
  smodes <- unique(smodes)
  if (length(smodes) != 1L) {
    stop("More than one data type across dataset: ",
         paste(smodes, collapse=','))
  }
  stopifnot(smodes == storage_mode)

  ## Should we enforce valid variable names in sample columns?
  ## biobroom says yes, but for now we just force that first character is a
  ## letter
  valid.colnames <- sapply(datasets, function(d) {
    all(substring(colnames(d), 1, 1) == make.names(substring(colnames(d), 1, 1)))
  })
  bad.ds <- which(!valid.colnames)
  if (length(bad.ds)) {
    stop("colnames(ds) should be valid variable names, offending entries:\n",
         paste(bad.ds, collapse=","))
  }

  ## Check that there are no duplicate rownames() in assays and that they all
  ## have measurements for the same features
  fids <- rownames(datasets[[1]])
  if (any(duplicated(fids))) stop("Duplicated rownames exist")
  features.consistent <- sapply(datasets, function(d) {
    nrow(d) == length(fids) && setequal(rownames(d), fids)
  })
  stopifnot(all(features.consistent))

  ## Check that we have feature_info entries for all features in assays
  ffi.equal <- setequal(facile_feature_info$feature_id, fids)
  if (!ffi.equal) {
    stop("facile_feature_info$feature_id is not seteqaul to datasets")
  }
  TRUE
}


#' Adds new set of assay data for all samples in a FacileDataSet
#'
#' Once a FacileDataSet has been created and initialized, either via a
#' low-level call to [initializeFacileDataSet()], or a call to
#' [as.FacileDataSet()] over a list of BiocAssayContainers, you can add more
#' assays (i.e. RNA-seq, microarray, etc.) to the FacileDataSet using this
#' function.
#'
#' Note that you cannot add assay data piecemeal. That is to say, you can not call
#' this function once to add copynumber data
#' (addFacileAssaySet(..., facile_assay_type = "cnv") to a subset of samples
#' and later call this function again to add copynumber to the rest of the
#' samples. The function will throw an error if
#' facile_assay_type %in% assay_names(x) == TRUE.
#'
#' @md
#' @importFrom rhdf5 h5createFile h5createDataset h5write
#' @export
#'
#' @param x The `FacileDataSet`
#' @param datasets list of `ExpressionSet`, `SummarizedExperiment`, or
#'   `DGEList`s that have the new assay data across all of the datasets in `x`.
#' @param facile_assay_name the name of the assay in the source dataset object
#' @param facile_assay_type string indicating the assay_type ('rnaseq',
#'   'affymetrix', etc.)
#' @param facile_feature_type a string indicating the universe the features in
#'   this can be anything, but certain identifiers are special like `"entrez"`,
#'   `"ensgid"`, `"enstid"`, etc.
#' @param facie_assay_description a string that allows the caller to provide
#'   a "freeform" description of the assay (platform, protocol, whatever).
#' @param facile_feature_info a `data.frame` with the required `feature_info`
#'   columns that describe the features in this assay. Please refer to the
#'   "Features" section of the `FacileDataSet` vignette for more complete
#'   description.
#' @param storage_mode either `"integer"` or `"numeric"`, maps to the
#'   `storage.mode` parameter in [rhdf5::h5createDataset()]
#' @param chunk_rows the first entry in the `chunk` parameter in
#'   [rhdf5::h5createDataset()] (`integer`)
#' @param chunk_cols the second entry in the `chunk` parameter in
#'   [rhdf5::h5createDataset()]. If this is `"ncol"`, it is set to the number
#'   of columns in each of the internal dataset matrices being added.
#' @param chunk_compression the `level` parameter in [rhdf5::h5createDataset()]
#' @param assay_name the assay name in the data containers provided in the
#'   `datasets` list.
#' @return a `tibble` subset of `facile_feature_info` that indicates the *new*
#'   features that were added to the internal `feature_info_tbl`.
#' @importFrom edgeR calcNormFactors
addFacileAssaySet <- function(x, datasets, facile_assay_name,
                              facile_assay_type,
                              facile_feature_type = NULL,
                              facile_assay_description = NULL,
                              facile_feature_info,
                              storage_mode = .storage.modes,
                              chunk_rows = 5000, chunk_cols = "ncol",
                              chunk_compression = 4,
                              assay_name = NULL, warn_existing = FALSE,
                              add_sample_covariates = TRUE,
                              covariate_def = NULL) {
  ## Parameter Checking --------------------------------------------------------
  stopifnot(is.FacileDataSet(x))
  assert_string(facile_assay_name)
  assert_string(facile_assay_type)
  if (facile_assay_name %in% assay_names(x)) {
    stop("`", facile_assay_name, "` assay already stored in FacileDataSet")
  }

  ainfo <- assay_info(x)
  if (facile_assay_type %in% ainfo$assay_type) {
    warning("assay exists of this type already: ", facile_assay_type,
            immediate.=TRUE)
  }
  facile_feature_type <- assert_string(facile_feature_type)
  if (is.null(facile_assay_description)) {
    facile_assay_description <- facile_assay_type
  }
  assert_string(facile_assay_description)
  storage_mode <- match.arg(storage_mode, .storage.modes)
  assert_valid_assay_datasets(datasets, facile_feature_info, storage_mode, assay_name)

  ## Ensure no redunancy in facile_feature_info
  nf <- facile_feature_info |> 
    distinct(feature_type, feature_id) |> 
    nrow()
  if (nrow(facile_feature_info) != nf) {
    stop("Redunant features in facile_feature_info")
  }

  ## Insert entry into assay_info table ----------------------------------------
  ai <- tibble(assay=facile_assay_name,
               assay_type=facile_assay_type,
               feature_type=facile_feature_type,
               description=facile_assay_description,
               nfeatures=nrow(datasets[[1]]),
               storage_mode=storage_mode) |>
    append_facile_table(x, "assay_info", warn_existing = warn_existing)

  ## Insert Feature Information into FacileDataSet -----------------------------
  ## Insert new features into global feature_info table
  features <- x |>
    append_facile_feature_info(facile_feature_info,
                               type = facile_feature_type,
                               warn_existing = warn_existing) |>
    select(feature_type, feature_id, added)

  ## Create entries in `assay_feature_info` table to track hdf5 indices for
  ## the features in this assay
  afi <- feature_info_tbl(x) |>
    semi_join(select(features, -added),
              by = c('feature_type', 'feature_id'),
              copy = TRUE, auto_index = TRUE) |>
    collect(n=Inf)
  afi <- transmute(afi, assay = facile_assay_name, feature_id,
                   hdf5_index = seq(nrow(afi))) |>
    append_facile_table(x, "assay_feature_info", warn_existing) |>
    arrange(hdf5_index)
  stopifnot(nf == nrow(afi), all(afi$hdf5_index == seq(nf)))

  ## Insert assay data and track sample info -----------------------------------
  aname <- paste0('assay/', facile_assay_name)
  stopifnot(h5createGroup(hdf5fn(x), aname))

  ## loop over datasets and insert one-by-one
  ## Mash it up all together (memory intensive!) to create norm.factor
  ## this also checks that rows are consistent across datasets (again) and
  ## returns them in the right order
  dats <- lapply(datasets, function(dat) {
    dat <- extract.assay(dat, assay_name)
    if (!setequal(rownames(dat), afi$feature_id)) {
      stop("Mismatch in rownames(", ds, ") and feature_id's in database")
    }
    dat[afi$feature_id,,drop=FALSE]
  })

  ## If rnaseq, do calcnormfactors
  if (facile_assay_type %in% c('rnaseq', 'isoseq', 'pseudobulk')) {
    ## Can create assay_sample_info_table
    cnts <- do.call(cbind, dats)
    normfactors <- calcNormFactors(cnts)
    asi <- tibble(
      assay=facile_assay_name,
      dataset=rep(names(datasets), sapply(dats, ncol)),
      sample_id=colnames(cnts),
      hdf5_index=lapply(dats, function(d) seq(ncol(d))) |> unlist(),
      libsize=colSums(cnts),
      normfactor=normfactors)
    rm(cnts)
  } else {
    asi <- lapply(names(dats), function(ds) {
      mtrx <- dats[[ds]]
      tibble(
        assay=facile_assay_name,
        dataset=ds,
        sample_id=colnames(mtrx),
        hdf5_index=seq(ncol(mtrx)),
        libsize=-1, normfactor=-1)
    }) |> bind_rows()
  }

  asi <- append_facile_table(asi, x, 'assay_sample_info', warn_existing)

  # "numeric" R storage mode is "double" storage mode in hdf5
  if (storage_mode == 'numeric') {
    storage_mode <- 'double'
  }

  assay.dat <- lapply(names(dats), function(ds) {
    ## This should be addFacileAssay(x, assay_name, dat, subset(asi, ...))
    dat <- dats[[ds]]
    if (is.character(chunk_cols) && chunk_cols == 'ncol') {
      xchunk.cols <- ncol(dat)
    } else if (is.integer(chunk_cols)) {
      xchunk.cols <- min(chunk_cols, ncol(dat))
    } else {
      stop("Unrecognized value for chunk.cols: ", chunk_cols)
    }
    xchunk.rows <- min(chunk_rows=5000, nrow(dat))
    chunk <- c(xchunk.rows, xchunk.cols)
    dname <- sprintf('%s/%s', aname, ds)
    if (chunk_compression == 0) {
      h5createDataset(hdf5fn(x), dname, dim(dat), storage.mode = storage_mode,
                      level = 0, filter = "NONE")
    } else {
      h5createDataset(hdf5fn(x), dname, dim(dat), storage.mode=storage_mode,
                      chunk = chunk, level = chunk_compression)
    }
    ## For some reason dispatching on this keeps messing up:
    ## Error in UseMethod("h5write") (from construction.R#296) :
    ##   no applicable method for 'h5write' applied to an object of class
    ##   "c('matrix', 'integer', 'numeric')"
    rhdf5::h5write(dat, file = hdf5fn(x), name = dname)
    tibble(assay=facile_assay_name, dataset=ds, sample_id=colnames(dat),
           hdf5_index=seq(ncol(dat)))
  })
  assay.dat <- bind_rows(assay.dat)

  ## Can check that asi == assay.dat
  chk <- inner_join(asi, assay.dat, by=c('assay', 'dataset', 'sample_id'))
  stopifnot(nrow(chk) == nrow(asi), all(chk$hdf5_index.x == chk$hdf5_index.y))

  # add samples to sample_info table if they're not there.
  samples <- asi |>
    mutate(parent_id=NA_character_) |>
    append_facile_table(x, 'sample_info', warn_existing)
  
  if (add_sample_covariates) {
    sample_covariates <- .prepare_sample_covariates(
      datasets, covariate_def = covariate_def)
    # add sample covariates to table
    sample.covs <- sample_covariates$eav |>
      mutate(date_entered = as.integer(Sys.time())) |>
      append_facile_table(x, "sample_covariate", 
                          warn_existing = warn_existing)
  }
  
  # TODO: Add sample covariates for samples that aren't yet loaded -------------
  invisible(list(samples=samples, assay_sample_info=asi))
}

#' Appends new features to `feature_info` table
#'
#' This function only adds features (feature_type, feature_id) that are not
#' in the `feature_info` table already
#'
#' @export
#' @param x The `FacileDataSet`
#' @param feature_info a table of new features that provides all columns
#'   in `feature_info_tbl(x)`
#' @param type A way to override (or set) the `feature_type` column of the
#'   `feature_info` table
#' @return invisible returns an annotated version of the `feature_info`
#'   table with an `$added` column with `TRUE/FALSE` values for the
#'   features that were new (and added) to the repository or `FALSE` to
#'   indicate that they were already in the database.
append_facile_feature_info <- function(x, feature_info,
                                       type = feature_info$feature_type,
                                       warn_existing = FALSE) {
  ## Argument Checking
  stopifnot(is.FacileDataSet(x))
  stopifnot(is.data.frame(feature_info))
  if (is.null(feature_info$feature_type)) {
    stopifnot(is.character(type), length(type) %in% c(1L, nrow(feature_info)))
    feature_info$feature_type <- type
  }
  ftypes <- unique(feature_info$feature_type)
  if (length(ftypes) > 1) {
    warning("Adding more than one feature_type to feature_info table",
            immediate.=TRUE)
  }

  # "source" may not be defined. If not, let's add it as "unspecified" and warn
  # the user
  if (is.null(feature_info[["source"]])) {
    warning('`"source"` column not specified in feature info, filling in with ',
            '"unspecified"', immediate. = TRUE)
    feature_info[["source"]] <- rep("unspecifed", nrow(feature_info))
  }
  
  added <- feature_info |>
    distinct(feature_type, feature_id, .keep_all = TRUE) |>
    append_facile_table(x, 'feature_info', warn_existing = warn_existing)
  invisible(added)
}
