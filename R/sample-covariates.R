##' Fetch rows from sample_covariate table for specified samples and covariates
##'
##' @export
##' @param db a \code{FacileDb} connection
##' @param samples a samples descriptor \code{tbl_*}
##' @param covariates character vector of covariate names
##' @param type only fetch covariates of a particular type?
##' @param do.collect for collection of result from database?
##' @return rows from the \code{sample_covaraite} table
fetch_sample_covariates <- function(db, samples=NULL, covariates=NULL) {
  stopifnot(is.FacileDb(db))
  dat <- sample_covariate_tbl(db)
  if (is.character(covariates)) {
    if (length(covariates) == 1L) {
      dat <- filter(dat, variable == covariates)
    } else if (length(covariates) > 1L) {
      dat <- filter(dat, variable %in% covariates)
    }
  }

  ## if (right.join) {
  ##   assert_sample_subset(samples)
  ##   out <- right_join(dat, samples, by=c('dataset', 'sample_id'))
  ## } else {
  ##   out <- filter_samples(dat, samples)
  ## }

  out <- filter_samples(dat, samples)
  set_fdb(out, db)
}

##' Appends covariate columns to a query result
##'
##' Note that this function will force the collection of \code{x}
##'
##' @export
##' @param x a "tbl"-like object that is a sample descriptor
with_sample_covariates <- function(x, covariates=NULL, db=fdb(x),
                                   cov.def=db[['cov.def']]) {
  stopifnot(is.FacileDb(db))
  stopifnot(is.character(covariates))

  samples <- assert_sample_subset(x) %>%
    select(dataset, sample_id) %>%
    collect(n=Inf) %>% ## can't call distinct on SQLite backend :-(
    distinct(.keep_all=TRUE)

  covs <- fetch_sample_covariates(db, samples, covariates) %>%
    spread_covariates(cov.def)

  collect(x, n=Inf) %>%
    left_join(covs, by=c('dataset', 'sample_id')) %>%
    set_fdb(db)
}

##' Spreads the covariates returned from database into wide data.frame
##'
##' Samples that did not have a value for a specific covariate are assigned to
##' have NA.
##'
##' @export
##' @param x output from \code{fetch_sample_covariates}
##' @return a wide \code{tbl_df}-like object
spread_covariates <- function(x, cov.def=NULL) {
  db <- NULL
  if (missing(cov.def) && is.null(cov.def) && is(x, 'tbl_sqlite')) {
    db <- fdb(x)
    if (is.FacileDb(db)) {
      cov.def <- db[['cov.def']]
    }
  }
  x <- assert_sample_covariates(x) %>%
    collect(n=Inf)

  ## Ensures we get a row for every sample in x, even if it is missing a value
  ## for the covariate
  dummy <- select(x, dataset, sample_id) %>%
    distinct(.keep_all=TRUE) %>%
    mutate(variable='.dummy.', value=NA)

  out <- bind_rows(x, dummy) %>%
    dcast(dataset + sample_id ~ variable, value.var='value') %>%
    mutate(.dummy.=NULL) %>%
    set_rownames(., paste(.$dataset, .$sample_id, sep='_'))

  if (!is.null(cov.def)) {
    for (cname in setdiff(colnames(out), c('dataset', 'sample_id'))) {
      casted <- cast_covariate(cname, out[[cname]], cov.def)
      if (is.data.frame(casted)) {
        out[[cname]] <- NULL
        out <- bind_cols(out, casted)
      } else {
        out[[cname]] <- casted
      }
      ## casting a survival covariate will return a two column thing with time
      ## and censoring information, so we need to account for that.
    }
  }

  set_fdb(out, db)
}

##' Casts the character values of the covariates to their defined types.
##'
##' For most things, a single value will be returned from each cast, but in the
##' case of "time_to_event" data, the value is expended to a two column
##' data.frame with a \code{tte_<covariate>} column for time to event, and an
##' \code{event_<covariate>} column to indicate event (1) or right censored (2).
##'
##' @importFrom yaml yaml.load_file
##' @export
##' @param covariate the name of the covariate
##' @param values the covariate values (which is a \code{character}) as it is
##'   pulled from the database.
##' @param cov.defs the path to the covariate definition yaml file
##'
##' @return values cast to appropriate type if a valid definition was found for
##'   \code{covariate}, otherwise values is returned "as is". Most of the time
##'   this is a single vector, but others it can be a data.frame (for
##'   \code{right_censored} data, for instance)
cast_covariate <- function(covariate, values, cov.def) {
  stopifnot(is.character(values))
  stopifnot(is.character(covariate) && length(covariate) == 1L)
  if (is.character(cov.def)) {
    cov.def <- yaml.load_file(assert_file(cov.def, 'r'))
  }
  stopifnot(is.list(cov.def))
  def <- cov.def[[covariate]]
  if (is.list(def)) {
    if (def$type == 'real') {
      values <- as.numeric(values)
    }
    if (def$type == 'right_censored') {
      values <- decode_right_censored(values, suffix=covariate)
    }
    if (is.character(def$levels)) {
      values <- factor(values, def$levels)
    }
  } else {
    if (covariate != '.dummy.') {
      warning("No covariate definition found for: ", covariate, immediate.=TRUE)
    }
  }

  values
}
