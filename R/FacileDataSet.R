#' Instantiates a FacileDataSet object from disk.
#'
#' The `FacileDataSet` is a reference data storage implementation that
#' implements the **FacileData Access API**. It facilitates the storage and
#' retrieval of large amounts of data by leveraging a SQLite database to store
#' sample- and feature-level metadata ("`pData`" and "`fData`"), and an HDF5
#' file to store all of the dense assay (matrix) data (gene counts, microarray
#' intensities, etc.) run over the samples.
#'
#' A `FacileDataSet` is materialized on disk by a well-structured directory,
#' which minimally includes the following items:
#'
#' 1. A `data.sqlite` SQLite database that stores feature and sample metadata
#' 2. A `data.h5` HDF5 file that stores a multitude of dense assay matrices that
#'    are generated from the assays performed on the samples in the
#'    `FacileDataSet`.
#' 3. A `meta.yaml` file tha contains informaiton about the `FacileDataSet`.
#'    To better understand the structure and contents of this file, you can
#'    refer to the following:
#'     a. The included `testdata/expected-meta.yaml` file, which is an
#'        exemplar file for [exampleFacileDataSet()].
#'     b. The help file provided by the [eav_metadata_create()] function, which
#'        describes in greater detail how we track a dataset's sample-level
#'        covariates (aka, "pData" in the bioconductor world).
#'    In the meantime, a short description of the entries found in the
#'    `meta.yaml` file is provded here:
#'     - `name`: the name of the dataset (ie. `"FacileTCGADataSet"`)
#'     - `organism`: `"Homo sapiens"`, `"Mus musculus"`, ec.
#'     - `default_assay`: the name of the assay to use by default if none is
#'       specified in calls to [fetch_assay_data()], [with_assay_data()], etc.
#'       (kind of like how `"exprs"` is the default assay used when working with
#'       a `Biobase::ExpressionSet`)
#'     - `datasets`: a section tha enumerates the datases included internally.
#'       The datasets are further enumerated.
#'     - `sample_covariates`: a section that enumerates the covariatets that
#'       are tracked over the samples inside the `FacileDataSet` (ie. a mapping
#'       of the `pData` for the samples). Reference [eav_metadata_create()]
#'       for more information.
#' 4. A `custom-annotation` directory, which stores custom `sample_covariate`
#'    (aka "pData") informaiton that analysts can identify and describe during
#'    the course of an analysis, or even add from external sources. Although
#'    this directory is required in the directory structure of a valid
#'    `FacileDataSet`, the `FacileDataSet()` constructor can be called with
#'    a custom `anno.dir` parameter so that custom annotations are stored
#'    elsewhere.
#'
#' @export
#' @importFrom RSQLite SQLite sqliteCopyDatabase
#' @importFrom DBI dbConnect dbDisconnect dbExecute
#' @family FacileDataSet
#'
#' @param path The path to the FacileData repository
#' @param data.fn A custom path to the database (probably don't mess with this)
#' @param sqlite.fn name of SQLite data file in FacileDataSet
#' @param hdf5.fn name of HDF5 data file in FacileDataSet
#' @param meta.fn name of metadata YAML data file in FacileDataSet
#' @param anno.dir A directory to house custom annotations/sample covariates
#' @param cache_size A custom paramter for the SQLite database
#' @param db.loc single character, location for the data
#' @param ... other args to pass down, not used at the moment
#' @param covdef.fn A custom path to the yaml file that has covariate mapping info
#' @return a `FacileDataSet` object
#' @examples
#' fn <- system.file("extdata", "exampleFacileDataSet", package = "FacileData")
#' fds <- FacileDataSet(fn)
FacileDataSet <- function(path, data.fn=file.path(path, 'data.sqlite'),
                          sqlite.fn=file.path(path, 'data.sqlite'),
                          hdf5.fn=file.path(path, 'data.h5'),
                          meta.fn=file.path(path, 'meta.yaml'),
                          anno.dir=file.path(path, 'custom-annotation'),
                          cache_size=80000,
                          db.loc=c('reference', 'temporary', 'memory'),
                          ...) {
  paths <- validate.facile.dirs(path, data.fn, sqlite.fn, hdf5.fn, meta.fn,
                                anno.dir)
  db.loc <- match.arg(db.loc)
  ## Update some parameters in the connection for increased speed
  ## http://stackoverflow.com/questions/1711631
  ##
  ## Explains some pragma setting:
  ##   http://www.codificar.com.br/blog/sqlite-optimization-faq/
  ##
  ## The following PRAGMA are of likely interest:
  ##   1. cache_size  https://www.sqlite.org/pragma.html#pragma_cache_size
  ##   2. page_size
  if (db.loc == 'temporary') {
    tmp.fn <- tempfile()
    file.copy(paths$sqlite.fn, tmp.fn)
    paths$sqlite.fn <- tmp.fn
  }

  con <- dbConnect(SQLite(), paths$sqlite.fn)

  if (db.loc == 'memory') {
    mcon <- dbConnect(RSQLite::SQLite(), ":memory:")
    RSQLite::sqliteCopyDatabase(con, mcon)
    dbDisconnect(con)
    con <- mcon
  }

  dbExecute(con, 'pragma temp_store=MEMORY;')
  dbExecute(con, sprintf('pragma cache_size=%d;', cache_size))

  if (!dir.exists(anno.dir)) {
    if (!dir.create(anno.dir, recursive =  TRUE))
      stop(sprintf("%s: failed to create annotation directory: ", FDS.name),
           anno.dir)
  }

  out <- list(con=con)
  out['parent.dir'] <- paths$path
  out['data.fn'] <- paths$sqlite.fn ## paths$data.fn
  out['sqlite.fn'] <- paths$sqlite.fn
  out['hdf5.fn'] <- paths$hdf5.fn
  out['anno.dir'] <- paths$anno.dir
  out['db.loc'] <- db.loc
  out[['cache']] <- new.env()

  ## meta information
  class(out) <- c("FacileDataSet", "FacileDataStore")

  mi <- meta_info(out)
  out['organism'] <- mi$organism

  if (is.null(mi$default_assay)) {
    mi$default_assay <- assay_names(out)[1L]
  }
  out['default_assay'] <- mi$default_assay

  class(out) <- c(mi$name, class(out))
  out
}

#' Class and validity checker for FacileDataSet
#'
#' @export
#'
#' @param x object to test
#' @return `TRUE`/`FALSE` indicating that `x` nominally "looks like" a
#'   `FacileDataSet`
is.FacileDataSet <- function(x) {
  test_facile_data_set(x)
}

#' Get location of the FacileDataSet database
#'
#' @export
#' @family FacileDataSet
#'
#' @param x FacileDataSet
#' @param mustWork boolean, if `TRUE` (default), throws an error if the sqlite
#'   file does not exist. When `FALSE`, this returns the "expected" path to the
#'   sqlite file for `x`
#' @return the filepath to the sqlite database
dbfn <- function(x, mustWork=TRUE) {
  base.fn <- 'data.sqlite'
  if (is.FacileDataSet(x)) {
    x <- x$parent.dir
  }
  assert_string(x)
  out <- file.path(x, base.fn)
  if (mustWork && !file.exists(out)) {
    stop("data.sqlite file not found: ", out)
  }
  out
}

#' Get location of the FacileDataSet HDF5 file
#'
#' @export
#' @family FacileDataSet
#'
#' @param x FacileDataSet
#' @param mustWork single logical
#' @return path to HDF5 file
hdf5fn <- function(x, mustWork=TRUE) {
  base.fn <- 'data.h5'
  if (is.FacileDataSet(x)) {
    x <- x$parent.dir
  }
  stopifnot(is.character(x), length(x) == 1L)
  out <- file.path(x, base.fn)
  if (mustWork && !file.exists(out)) {
    stop("data.h5 file not found: ", out)
  }
  out
}

#' Path to the meta information YAML file
#'
#' @export
#' @rdname meta-info
#' @family FacileDataSet
#'
#' @param x A `FacileDataSet`
meta_file <- function(x) {
  assert_facile_data_set(x)
  fn <- assert_file(file.path(x$parent.dir, 'meta.yaml'), 'r')
  fn
}

#' Retrieves the meta information for a FacileDataSet
#'
#' Lots of useful information is stored in a `FacileDataSet`'s `meta.yaml` file.
#' This function returns all of that in a list-of-lists
#'
#' @export
#' @rdname meta-info
#' @param fn The path to the `meta.yaml` file.
#' @return The `meta.yaml` file parsed into a list-of-lists representation
meta_info <- function(x, fn = meta_file(x)) {
  assert_facile_data_set(x)
  out <- assert_valid_meta_file(fn, as.list = TRUE)
  out
}

#' Retrieves the organism the data is defined over
#'
#' A FacileDataStore is only expected to hold data for one organism.
#'
#' @export
#' @family API
#' @return `"Homo sapiens`", `"Mus musculus"`, etc.
organism.FacileDataSet <- function(x) {
  assert_facile_data_set(x)
  x$organism
}

#' @export
#' @rdname meta-info
default_assay.FacileDataSet <- function(x) {
  if (is.null(x$default_assay)) {
    out <- assay_names(x, default_first=FALSE)[1L]
    if (is.na(out)) out <- NULL
  } else {
    out <- x$default_assay
  }
  out
}

#' Retrieves URL and description of datasets in a FacileDataSet
#'
#' A `FacileDataSet` can contain assay data from different "datasets" (such
#' as different cancer indications from the TCGA). This functions returns
#' description and URL information that describes these datasets in more detail,
#' which is specified in the FacileDataSets `meta.yaml` file.
#'
#' @rdname meta-info
#' @export
#' @param as.list boolean, if `FALSE` (default) returns a list, otherwise
#'   summarizes results into a tibble.
#' @return meta information about the datasets in `x` as a `list` or `tibble`
dataset_definitions <- function(x, as.list=TRUE) {
  defs <- meta_info(x)$datasets
  if (!as.list) {
    defs <- lapply(names(defs), function(ds) {
      i <- defs[[ds]]
      tibble(dataset=ds, url=i$url, description=i$description)
    })
    defs <- bind_rows(defs)
  }
  defs
}

#' Get description of sample metadata columns
#'
#' Descriptions of the sample covariates can be specified in a FacileDataSet's
#' `meta.yaml` file. This function returns those.
#'
#' @export
#' @importFrom yaml yaml.load_file
#' @param x FacileDataSet
#' @param as.list single logical, return tibble or list
#' @return meta information about the sample covariates in `x`
covariate_definitions <- function(x, as.list=TRUE) {
  out <- meta_info(x)$sample_covariates
  if (!as.list) {
    out <- lapply(names(out), function(name) {
      i <- out[[name]]
      lvls <- unique(i$levels)
      is.factor <- !is.null(lvls)
      lbl <- if (is.null(i$label)) name else i$label
      tibble(variable=name, type=i$type, class=i$class, label=i$label,
             is_factor=is.factor, levels=list(lvls), description=i$description)
    })
    out <- bind_rows(out)
  }
  class(out) <- c('CovariateDefinitions', class(out))
  set_fds(out, x)
}

#' Retrieves the sample identifiers for all samples in a FacileDataSet.
#'
#' Sample identifiers are provided as `dataset,sample_id tuples`.
#'
#' @export
#' @family API
#'
#' @param x a `FacileDataSet`
#' @return tibble of sample attributes
samples.FacileDataSet <- function(x, ...) {
  assert_facile_data_set(x)
  sample_info_tbl(x) %>%
    select(dataset, sample_id) %>%
    as_facile_frame(x, ...)
}

#' @export
#' @rdname facet_frame
facet_frame.FacileDataSet <- function(x, name = "default", ...) {
  samples(x) %>%
    mutate(facet = dataset) %>%
    select(facet, dataset, sample_id) %>%
    as_facile_frame(x)
}

#' @noRd
#'
#' @export
#' @importFrom DBI dbGetQuery
#' @param x FacileDataSet
#' @param ... additional args (ignored)
print.FacileDataSet <- function(x, ...) {
  ns <- dbGetQuery(x$con, "SELECT COUNT(*) FROM sample_info;")
  nds <- dbGetQuery(x$con, "SELECT COUNT (DISTINCT dataset) FROM sample_info;")
  out <- paste0(
    "FacileDataSet",
    if (class(x)[1] != "FacileDataSet") sprintf(" (%s)\n", class(x)[1L]) else "\n",
    "  Directory: ", x$parent.dir, "\n",
    sprintf("  %d samples over %d datasets\n", ns[1,1], nds[1,1]))
  cat(out)
  invisible()
}
