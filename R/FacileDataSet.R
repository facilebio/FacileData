##' Connect to a FacileDataSet repository.
##'
##' @export
##' @importFrom RSQLite dbConnect SQLite dbExecute
##' @param path The path to the FacileData repository
##' @param data.fn A custom path to the database (probably don't mess with this)
##' @param covdef.fn A custom path to the yaml file that has covariate mapping info
##' @param anno.dir A directory to house custom annotations/sample covariates
##' @param cache_size A custom paramter for the SQLite database
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
  ## Explains some pragme setting:
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
    RSQLite::dbDisconnect(con)
    con <- mcon
  }

  dbExecute(con, 'pragma temp_store=MEMORY;')
  dbExecute(con, sprintf('pragma cache_size=%d;', cache_size))

  out <- list(con=con)
  out['parent.dir'] <- paths$path
  out['data.fn'] <- paths$sqlite.fn ## paths$data.fn
  out['sqlite.fn'] <- paths$sqlite.fn
  out['hdf5.fn'] <- paths$hdf5.fn
  out['anno.dir'] <- paths$anno.dir
  out['db.loc'] <- db.loc
  out[['cache']] <- list()

  ## meta information
  class(out) <- 'FacileDataSet'

  mi <- meta_info(out)
  out['organism'] <- mi$organism

  if (is.null(mi$default_assay)) {
    mi$default_assay <- assay_names(out)[1L]
  }
  out['default_assay'] <- mi$default_assay

  class(out) <- c(mi$name, class(out))

  out
}

##' @export
is.FacileDataSet <- function(x) {
  ## `is` and `validate` funcitonality confused in here.
  is(x, 'FacileDataSet') &&
    'con' %in% names(x) && is(x$con, 'DBIObject') &&
    'anno.dir' %in% names(x) && dir.exists(x$anno.dir) &&
    'hdf5.fn' %in% names(x) && file.exists(x$hdf5.fn)
}

dbfn <- function(x, mustWork=TRUE) {
  base.fn <- 'data.sqlite'
  if (is.FacileDataSet(x)) {
    x <- x$parent.dir
  }
  stopifnot(is.character(x), length(x) == 1L)
  out <- file.path(x, base.fn)
  if (mustWork && !file.exists(out)) {
    stop("data.sqlite file not found: ", out)
  }
  out
}

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

##' Path to the meta information file
##'
##' @rdname meta-info
##' @export
##' @param x \code{FacileDataSet}
meta_file <- function(x) {
  stopifnot(is.FacileDataSet(x))
  fn <- assert_file(file.path(x$parent.dir, 'meta.yaml'), 'r')
  fn
}

##' Get meta information for dataset
##'
##' @rdname meta-info
##' @export
##' @param x \code{FacileDataSet}
meta_info <- function(x) {
  out <- yaml.load_file(meta_file(x))
  out
}

##' @rdname meta-info
##' @export
##' @param x \code{FacileDataSet}
organism <- function(x) {
  stopifnot(is.FacileDataSet(x))
  x$organism
}

##' @rdname meta-info
##' @export
##' @param x \code{FacileDataSet}
default_assay <- function(x) {
  stopifnot(is.FacileDataSet(x))
  if (is.null(x$default_assay)) {
    out <- assay_names(x, default_first=FALSE)[1L]
    if (is.na(out)) out <- NULL
  } else {
    out <- x$default_assay
  }
  out
}

##' @rdname meta-info
##' @export
dataset_definitions <- function(x, as.list=TRUE) {
  defs <- meta_info(x)$datasets
  if (!as.list) {
    defs <- lapply(names(defs), function(ds) {
      i <- defs[[ds]]
      tibble(dataset=ds, url=i$url, description=i$description)
    }) %>% bind_rows
  }
  defs
}

##' @export
##' @importFrom yaml yaml.load_file
covariate_definitions <- function(x, as.list=TRUE) {
  out <- meta_info(x)$sample_covariates
  if (!as.list) {
    out <- lapply(names(out), function(name) {
      i <- out[[name]]
      lvls <- i$levels
      is.factor <- !is.null(lvls)
      lbl <- if (is.null(i$label)) name else i$label
      tibble(variable=name, type=i$type, class=i$class, label=i$label,
             is_factor=is.factor, levels=list(lvls), description=i$description)
    }) %>% bind_rows
  }
  class(out) <- c('CovariateDefinitions', class(out))
  out %>% set_fds(x)
}

##' @export
samples <- function(x) {
  stopifnot(is.FacileDataSet(x))
  sample_info_tbl(x) %>%
    select(dataset, sample_id) %>%
    set_fds(x)
}

##' A dataset-specific method to overridden that returns default sample grouping
##'
##' @export
##' @rdname facet_frame
facet_frame <- function(x, ...) {
  UseMethod("facet_frame")
}

##' @export
##' @rdname facet_frame
facet_frame.default <- function(x, ...) {
  tibble(facet=character(), dataset=facet, sample_id=facet)
}

##' @export
print.FacileDataSet <- function(x, ...) {
  ns <- RSQLite::dbGetQuery(x$con, "SELECT COUNT(*) FROM sample_info;")
  nds <- RSQLite::dbGetQuery(x$con, "SELECT COUNT (DISTINCT dataset) FROM sample_info;")
  out <- paste0(
    "FacileDataSet",
    if (class(x)[1] != "FacileDataSet") sprintf(" (%s)\n", class(x)[1L]) else "\n",
    "  Directory: ", x$parent.dir, "\n",
    sprintf("  %d samples over %d datasets\n", ns[1,1], nds[1,1]))
  cat(out)
  invisible()
}

