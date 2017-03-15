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

  # out['parent.dir'] <- list(NULL)
  # out['data.fn'] <- list(NULL)
  # out['cov.def'] <- list(NULL)
  # out['anno.dir'] <- list(NULL)
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

##' @export
covariate_definition_file <- function(x) {
  stopifnot(is.FacileDataSet(x))
  x[['cov.def']]
}

##' @export
##' @importFrom yaml yaml.load_file
covariate_definitions <- function(x) {
  stopifnot(is.FacileDataSet(x))
  cov.def <- covariate_definition_file(x)
  assert_file(cov.def, 'r')
  out <- yaml.load_file(cov.def)
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

