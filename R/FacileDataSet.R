##' Connect to a FacileDataSet repository.
##'
##' @export
##' @param path The path to the FacileData repository
##' @param data.fn A custom path to the database (probably don't mess with this)
##' @param covdef.fn A custom path to the yaml file that has covariate mapping info
##' @param anno.dir A directory to house custom annotations/sample covariates
##' @param cache_size A custom paramter for the SQLite database
FacileDataSet <- function(path, data.fn=file.path(path, paste0('data.', db.type)),
                          covdef.fn=file.path(path, 'sample-covariate-info.yaml'),
                          anno.dir=file.path(path, 'custom-annotation'),
                          cache_size=80000, db.type=c('sqlite', 'monetdblite')) {
  db.type <- match.arg(db.type)
  paths <- validate.facile.dirs(path, data.fn, covdef.fn, anno.dir, db.type)

  if (db.type == 'sqlite') {
    ## Update some parameters in the connection for increased speed
    ## http://stackoverflow.com/questions/1711631
    ##
    ## Explains some pragme setting:
    ##   http://www.codificar.com.br/blog/sqlite-optimization-faq/
    ##
    ## The following PRAGMA are of likely interest:
    ##   1. cache_size  https://www.sqlite.org/pragma.html#pragma_cache_size
    ##   2. page_size
    out <- src_sqlite(paths$data.fn)
    dbSendQuery(out$con, 'pragma temp_store=MEMORY;')
    dbSendQuery(out$con, sprintf('pragma cache_size=%d;', cache_size))
  } else if (db.type == 'monetdblite') {
    if (!require('MonetDBLite')) {
      stop("MonetDBLite required to access MonetDBLite database")
    }
    out <- src_monetdblite(paths$data.fn)
  }

  # out['parent.dir'] <- list(NULL)
  # out['data.fn'] <- list(NULL)
  # out['cov.def'] <- list(NULL)
  # out['anno.dir'] <- list(NULL)
  out['parent.dir'] <- paths$path
  out['data.fn'] <- paths$data.fn
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

##' @export
expression_tbl <- function(x) {
  stopifnot(is.FacileDataSet(x))
  tbl(x, 'expression') %>% set_fds(x)
}

##' @export
sample_stats_tbl <- function(x) {
  stopifnot(is.FacileDataSet(x))
  tbl(x, 'sample_stats') %>% set_fds(x)
}

##' @export
sample_covariate_tbl <- function(x) {
  stopifnot(is.FacileDataSet(x))
  tbl(x, 'sample_covariate') %>% set_fds(x)
}

##' @export
gene_info_tbl <- function(x) {
  stopifnot(is.FacileDataSet(x))
  tbl(x, 'gene_info') %>% set_fds(x)
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

##' Get/set db
##'
##' @rdname getsetdb
##' @export
##' @param x the object
##' @param db The \code{FacileDb} object
fds <- function(x) {
  out <- attr(x, 'fds')
  if (is.null(out)) {
    warning("No FacileDataSet found in x (", class(x)[1L], ")", immediate.=TRUE)
  }
  out
}

##' @rdname getsetdb
##' @export
"fds<-" <- function(x, value) {
  UseMethod("fds<-", x)
}

##' @rdname getsetdb
##' @export
"fds<-.tbl" <- function(x, value) {
  attr(x, 'fds') <- value
  x
}

##' @rdname getsetdb
##' @export
"fds<-.data.frame" <- function(x, value) {
  attr(x, 'fds') <- value
  x
}

"fds<-.default" <- function(x, value) {
  attr(x, 'fds') <- value
  x
}

##' @rdname getsetdb
##' @export
set_fds <- function(x, value) {
  attr(x, 'fds') <- value
  x
}

## Unexported utility functions ================================================

validate.facile.dirs <- function(path, data.fn, covdef.fn, anno.dir,
                                 db.type=c('sqlite', 'monetdblite')) {
  if (!dir.exists(path)) {
    stop("Top level FacileData directory does not exist: ", path)
  }
  path <- normalizePath(path)
  if (!file.exists(data.fn)) {
    stop("Data file (database) does not exists", data.fn)
  } else {
    data.fn <- normalizePath(data.fn)
    if (dirname(data.fn) != path) {
      warning("Data file is not under parent directory", immediate.=TRUE)
    }
  }
  if (!file.exists(covdef.fn)) {
    stop("Covariate information yaml file does not exist: ", covdef.fn)
  } else {
    covdef.fn <- normalizePath(covdef.fn)
    if (dirname(covdef.fn) != path) {
      warning("Covariate file is not under parent directory", immediate.=TRUE)
    }
  }
  if (!dir.exists(anno.dir)) {
    stop("Directory for custom annotations does not exist: ", anno.dir)
  } else {
    anno.dir <- normalizePath(anno.dir)
    if (dirname(anno.dir) != path) {
      warning("Custom annotation directory not under parent directory.",
              immediate.=TRUE)
    }
  }
  db.type <- match.arg(db.type)

  list(path=path, data.fn=data.fn, covdef.fn=covdef.fn, anno.dir=anno.dir)
}
