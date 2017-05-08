##' Connect to a FacileDataSet repository.
##'
##' @export
##' @param path The path to the FacileData repository
##' @param data.fn A custom path to the database (probably don't mess with this)
##' @param covdef.fn A custom path to the yaml file that has covariate mapping info
##' @param anno.dir A directory to house custom annotations/sample covariates
##' @param cache_size A custom paramter for the SQLite database
FacileDataSet <- function(path, data.fn=file.path(path, 'data.sqlite'),
                          sqlite.fn=file.path(path, 'data.sqlite'),
                          hdf5.fn=file.path(path, 'data.h5'),
                          covdef.fn=file.path(path, 'sample-covariate-info.yaml'),
                          anno.dir=file.path(path, 'custom-annotation'),
                          cache_size=80000,
                          db.loc=c('reference', 'temporary', 'memory')) {
  paths <- validate.facile.dirs(path, data.fn, sqlite.fn, hdf5.fn, covdef.fn,
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
    out <- src_sqlite(paths$sqlite.fn)
  } else {
    out <- src_sqlite(paths$sqlite.fn)
    if (db.loc == 'memory') {
      mcon <- dbConnect(RSQLite::SQLite(), ":memory:")
      RSQLite::sqliteCopyDatabase(out$con, mcon)
      RSQLite::dbDisconnect(out$con)
      out$con <- mcon
    }
  }

  dbGetQuery(out$con, 'pragma temp_store=MEMORY;')
  dbGetQuery(out$con, sprintf('pragma cache_size=%d;', cache_size))

  out['parent.dir'] <- paths$path
  out['data.fn'] <- paths$sqlite.fn ## paths$data.fn
  out['sqlite.fn'] <- paths$sqlite.fn
  out['hdf5.fn'] <- paths$hdf5.fn
  out['cov.def'] <- paths$covdef.fn
  out['anno.dir'] <- paths$anno.dir

  class(out) <- c('FacileDataSet', class(out))
  out
}

##' @export
is.FacileDataSet <- function(x) {
  is(x, 'FacileDataSet') &&
    'con' %in% names(x) && is(x$con, 'DBIObject') &&
    'cov.def' %in% names(x) && file.exists(x$cov.def) &&
    'anno.dir' %in% names(x) && dir.exists(x$anno.dir)
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
    collect(n=Inf)
}
