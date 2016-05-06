##' Fetch rows from sample_covariate table for specified samples and covariates
##'
##' @export
##' @param db a \code{FacileDb} connection
##' @param samples a samples descriptor \code{tbl_*}
##' @param covariates character vector of covariate names
##' @param do.collect for collection of result from database?
##' @return rows from the \code{sample_covaraite} table
fetch_sample_covariates <- function(db, samples=NULL, covariates=NULL,
                                    do.collect=FALSE) {
  stopifnot(is.FacileDb(db))
  dat <- sample_covariate_tbl(db)
  if (is.character(covariates)) {
    if (length(covariates) == 1L) {
      dat <- filter(dat, variable == covariates)
    } else if (length(covariates) > 1L) {
      dat <- filter(dat, variable %in% covariates)
    }
  }
  out <- filter_samples(dat, samples)
  if (do.collect) {
    out <- collect(out)
    db <- NULL
  }
  attr(out, 'db') <- db
  out
}

##' Appends covariate columns to a query result
##'
##' Note that this function will force the collection of \code{x}
##'
##' @export
with_sample_covariates <- function(x, covariates=NULL, db=attr(x, 'db'),
                                   cov.def=db[['cov.def']]) {
  stopifnot(is.FacileDb(db))
  stopifnot(is.character(covariates))

  samples <- assert_sample_subset(x) %>%
    select(dataset, sample_id) %>%
    distinct

  covs <- fetch_sample_covariates(db, samples, covariates) %>%
    spread_covariates(cov.def)

  inner_join(collect(x), covs, by=c('dataset', 'sample_id'))
}

##' Spreads the covariates returned from database into wide data.frame
##'
##' Samples that did not have a value for a specific covariate are assigned to
##' have NA.
##'
##' @param x output from \code{fetch_sample_covariates}
##' @return a wide \code{tbl_df}-like object
spread_covariates <- function(x, cov.def=NULL) {
  if (missing(cov.def) && is.null(cov.def) && is(x, 'tbl_sqlite')) {
    db <- attr(x, 'db')
    if (is.FacileDb(db)) {
      cov.def <- db[['cov.def']]
    }
  }
  x <- assert_sample_covariates(x) %>%
    collect
  dummy <- select(x, dataset, sample_id) %>%
    distinct %>%
    mutate(variable='.dummy.', value=NA)
  out <- bind_rows(x, dummy) %>%
    dcast(dataset + sample_id ~ variable, value.var='value') %>%
    set_rownames(., paste(.$dataset, .$sample_id, sep='_'))
  if (!is.null(cov.def)) {
    for (cname in setdiff(colnames(out), c('dataset', 'sample_id'))) {
      out[[cname]] <- cast_covariate(cname, out[[cname]], cov.def)
    }
  }
  transform(out, .dummy.=NULL)
}

##' Casts the character values of the covariates to their defined types
##'
##' @importFrom yaml yaml.load_file
##' @export
##' @param covariate the name of the covariate
##' @param values the covariate values (likely as a character)
##' @param cov.defs the path to the covariate definition yaml file
##' @return values cast to appropriate type if a valid definition was found for
##'   \code{covariate}, otherwise values is returned "as is"
cast_covariate <- function(covariate, values, cov.def) {
  if (is.character(cov.def)) {
    cov.def <- yaml.load_file(assert_file(cov.def, 'r'))
  }
  stopifnot(is.list(cov.def))
  def <- cov.def[[covariate]]
  if (is.list(def)) {
    if (def$type == 'real') {
      values <- as.numeric(values)
    }
    if (is.character(def$levels)) {
      values <- factor(values, def$levels)
    }
  }
  values
}
