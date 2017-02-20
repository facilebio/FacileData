## The code in here builds FacileDataSets from Bioconductor like objects

##' Create an empty FacileDataSet
##'
##' @importFrom rhdf5 h5createFile h5createGroup H5close
##' @export
##'
##' @param path the directory to create which will house the
##'   \code{FacileDataSet}
##' @param covariate_definition the path to the covariate definition file
##' @param page_size,cache_size \code{pragma} values to setup the backend SQLite
##'   database
##' @return inivisibly returns the \code{FaclieDataSet} you just made
initializeFacileDataSet <- function(path, covariate_definition,
                                    page_size=2**12, cache_size=2e5) {
  path <- normalizePath(path, mustWork=FALSE)
  if (file.exists(path) || dir.exists(path)) {
    stop("Collision with file: ", path)
  }
  parent.dir <- dirname(path)
  assert_directory_exists(parent.dir)
  assert_access(parent.dir, 'w') ## check parent directory is writeable

  dir.create(path)
  dir.create(file.path(path, 'custom-annotation'))
  file.copy(covariate_definition, file.path(path, 'sample-covariate-info.yaml'))

  ## Create sqlite db
  db.fn <- file.path(path, 'data.sqlite')
  con <- dbConnect(SQLite(), db.fn)
  on.exit(dbDisconnect(con))
  dbGetQuery(con, 'pragma temp_store=MEMORY;')
  dbGetQuery(con, sprintf('pragma page_size=%d', page_size))
  dbGetQuery(con, sprintf('pragma cache_size=%d;', cache_size))
  sql.fn <- system.file('extdata', 'init', 'faciledataset.sql',
                        package='FacileDataSet')
  db.sql <- sqlFromFile(sql.fn)
  dbGetQueries(con, db.sql)

  ## Create empty HDF5 file
  hd5.fn <- file.path(path, 'data.h5')
  h5createFile(hd5.fn)
  h5createGroup(hd5.fn, 'assay')
  H5close()

  invisible(FacileDataSet(path))
}

.feature.types <- c('entrez', 'ensgid', 'enstid')
.assay.types <- c('rnaseq', 'isoseq', 'nanostring', 'fluidigm')
.storage.modes <- c('integer', 'numeric')

supported.assay.container <- function(x) {
  ## Currently only support DGEList
  is(x, 'DGEList') || is(x, 'eSet') || is(x, 'SummarizedExperiment')
}

extract.assay <- function(x, assay_name=NULL) {
  ## OO? We don't need no stinking OO
  if (is(x, 'DGEList')) {
    out <- x$counts
  } else if (is(x, 'eSet')) {
    ns <- loadNamespace("Biobase")
    if (is.null(assay_name)) assay_name <- assayDataElementNames(x)[1]
    out <- ns$assayDataElement(x, assay_name)
  } else if (is(x, 'SummarizedExperiment')) {
    ns <- loadNamespace("SummarizedExperiment")
    out <- if (is.null(assay_name)) assays(x)[[1L]] else assay(x, assay_name)
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
    bad <- sapply(datasets[idxs], function(d) class(d)[1L]) %>% unique
    stop("Unsupported assay container(s):", paste(bad, collapse=","))
  }
  dnames <- names(datasets)
  if (is.null(dnames) ||
      all(dnames != make.names(dnames)) ||
      length(unique(dnames)) != length(dnames)) {
    stop("names(datasets) must be valid varnames and all unique")
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
    nrow(d) == length(fids) & setequal(rownames(d), fids)
  })
  stopifnot(all(features.consistent))

  ## Check that we have feature_info entries for all features in assays
  ffi.equal <- setequal(facile_feature_info$feature_id, fids)
  if (!ffi.equal) {
    stop("facile_feature_info$feature_id is not seteqaul to datasets")
  }

  TRUE
}


##' Adds a complete set of assay data for all samples across datasets in the
##' FacileDataSet
##'
##' This needs to be busted up into functions Minimally loop over datasets to
##' addFacileAssay
##'
##' @importFrom rhdf5 h5createFile h5createDataset h5write
##' @export
##'
##' @param x The \code{FacileDataSeta}
##' @param dat list of ExpressionSet, SummarizedExperiment, or DGELists that
##'   have the assay data for the given assay across all of our datasets
##' @param assay_name the name of the assay in the source dataset object
##' @param facile_assay_name the name of the assay to store in the FacileDataSet
##' @param facile_assay_type string indicating the assay_type
addFacileAssaySet <- function(x, datasets, facile_assay_name,
                              facile_assay_type=.assay.types,
                              facile_feature_type=.feature.types,
                              facile_assay_description=NULL,
                              facile_feature_info,
                              storage_mode=.storage.modes,
                              chunk_rows=5000, chunk_cols='ncol',
                              chunk_compression=4,
                              assay_name=NULL) {
  ## Parameter Checking --------------------------------------------------------
  stopifnot(is.FacileDataSet(x))
  assert_string(facile_assay_name)
  if (facile_assay_name %in% assay_names(x)) {
    stop("`", facile_assay_name, "` assay already stored in FacileDataSet")
  }
  facile_assay_type <- match.arg(facile_assay_type, .assay.types)
  if (facile_assay_type %in% assay_types(x)) {
    warning("assay exists of this type already: ", facile_assay_type,
            immediate.=TRUE)
  }
  facile_feature_type <- match.arg(facile_feature_type, .feature.types)
  if (is.null(facile_assay_description)) {
    facile_assay_description <- facile_assay_type
  }
  assert_string(facile_assay_description)
  storage_mode <- match.arg(storage_mode, .storage.modes)
  assert_valid_assay_datasets(datasets, facile_feature_info, storage_mode)

  ## Insert entry into assay_info table ----------------------------------------
  ai <- tibble(assay=facile_assay_name,
               assay_type=facile_assay_type,
               feature_type=facile_feature_type,
               description=facile_assay_description,
               nfeatures=nrow(datasets[[1]]),
               storage_mode=storage_mode) %>%
    append_facile_table(x, 'assay_info')

  ## Insert Feature Information into FacileDataSet -----------------------------
  ## Insert new features into global feature_info table
  features <- append_facile_feature_info(x, facile_feature_info,
                                         type=facile_feature_type)

  ## Create entries in `assay_feature_info` table to track hdf5 indices for
  ## the features in this assay
  afi <- features %>%
    transmute(assay=facile_assay_name, feature_id, hdf5_index=seq(nrow(.))) %>%
    append_facile_table(x, 'assay_feature_info')

  stopifnot(nrow(features) == nrow(afi))

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

  ## Can create assay_sample_info_table
  y <- edgeR::DGEList(do.call(cbind, dats)) %>% calcNormFactors
  asi <- tibble(
    assay=facile_assay_name,
    dataset=rep(names(datasets), sapply(dats, ncol)),
    sample_id=colnames(y),
    hdf5_index=lapply(dats, function(d) seq(ncol(d))) %>% unlist,
    libsize=y$samples$lib.size,
    normfactor=y$samples$norm.factors) %>%
    append_facile_table(x, 'assay_sample_info')
  rm(y)
  gc()

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
    h5createDataset(hdf5fn(x), dname, dim(dat), storage.mode=storage_mode,
                    chunk=chunk, level=chunk_compression)
    h5write(dat, file=hdf5fn(x), dname)
    tibble(assay=facile_assay_name, dataset=ds, sample_id=colnames(dat),
           hdf5_index=seq(ncol(dat)))
  })
  assay.dat <- bind_rows(assay.dat)

  ## Can check that asi == assay.dat
  chk <- inner_join(asi, assay.dat, by=c('assay', 'dataset', 'sample_id'))
  stopifnot(nrow(chk) == nrow(asi), all(chk$hdf5_index.x == chk$hdf5_index.y))

  ## add samples to sample_info table if they're not there.
  samples <- asi %>%
    mutate(parent_id=NA_character_) %>%
    append_facile_table(x, 'sample_info')
  invisible(samples)
}

##' Appends new features to \code{feature_info} table
##'
##' This function only adds features (feature_type, feature_id) that are not
##' in the \code{feature_info} table already
##'
##' @export
##' @param x The \code{FacileDataSet}
##' @param feature_info a table of new features that provides all columns
##'   in \code{feature_info_tbl(x)}
##' @param type A way to override (or set) the \code{feature_type} column of the
##'   \code{feature_info} table
##' @return invisible returns an annotated version of the \code{feature_info}
##'   table with an \code{$added} column with \code{TRUE/FALSE} values for the
##'   features that were new (and added) to the repository or \code{FALSE} to
##'   indicate that they were already in the database.
append_facile_feature_info <- function(x, feature_info,
                                       type=feature_info$feature_type) {
  ## Argument Checking
  stopifnot(is.FacileDataSet(x))
  stopifnot(is.data.frame(feature_info))
  if (is.null(feature_info$feature_type)) {
    stopifnot(is.character(type), length(type) %in% c(1L, nrow(feature_info)))
    feature_info$feature_type <- type
  }
  ftypes <- unique(facile_feature_info$feature_type)
  if (length(ftypes) > 1) {
    warning("Adding more than one feature_type to feature_info table",
            immediate.=TRUE)
  }
  stopifnot(all(ftypes %in% .feature.types))
  added <- facile_feature_info %>%
    distinct(feature_type, feature_id, .keep_all=TRUE) %>%
    append_facile_table(x, 'feature_info')
  invisible(added)
}
