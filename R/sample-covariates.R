##' Fetch rows from sample_covariate table for specified samples and covariates
##'
##' @export
##' @param db a \code{FacileDataSet} connection
##' @param samples a samples descriptor \code{tbl_*}
##' @param covariates character vector of covariate names
##' @param custom_key The key to use to fetch more custom annotations over
##'   the given samples
##' @return rows from the \code{sample_covaraite} table
fetch_sample_covariates <- function(x, samples=NULL, covariates=NULL,
                                    custom_key=Sys.getenv("USER")) {
  stopifnot(is.FacileDataSet(x))
  dat <- sample_covariate_tbl(x)
  if (is.character(covariates)) {
    if (length(covariates) == 1L) {
      dat <- filter(dat, variable == covariates)
    } else if (length(covariates) > 1L) {
      dat <- filter(dat, variable %in% covariates)
    }
  }

  ## If the samples descriptor is defined over the sample_covariate table,
  ## this thing explodes (inner joining within itself, I guess). We defensively
  ## copy the sample descriptor, but in future maybe better to test if the
  ## dat and samples sqlite tables are pointing to the same thing
  if (!is.null(samples)) {
    assert_sample_subset(samples)
    samples <- collect(samples, n=Inf) %>% distinct(dataset, sample_id)
  }
  out <- filter_samples(dat, samples)

  if (!is.null(custom_key)) {
    custom <- fetch_custom_sample_covariates(x, samples,
                                             covariates=covariates,
                                             custom_key)
    out <- bind_rows(collect(out, n=Inf), custom)
  }

  set_fds(out, x)
}

##' Fetches custom (user) annotations for a given user prefix
##'
##' @export
##' @importFrom jsonlite stream_in
##' @param fds The \code{FacileDataSet}
##' @param samples the facile sample descriptor
##' @param custom_key The key to use for the custom annotation
##' @return covariate tbl
fetch_custom_sample_covariates <- function(x, samples=NULL, covariates=NULL,
                                           custom_key=Sys.getenv("USER"),
                                           file.prefix="facile") {
  stopifnot(is.FacileDataSet(x))
  out.cols <- colnames(sample_covariate_tbl(x))

  fpat <- paste0('^', file.prefix, '_', custom_key, "_*")
  annot.files <- list.files(path=x$anno.dir, pattern=fpat, full.names=TRUE)

  if (length(annot.files)) {
    annos <- lapply(annot.files, function(fn) stream_in(file(fn), verbose=FALSE))
    out <- bind_rows(annos) %>%
      mutate(class='user_annotation') %>%
      select_(.dots=out.cols) %>%
      set_fds(x) %>%
      filter_samples(samples)
    ## We weren't saving the type == 'categorical' column earlier. So if this
    ## column is.na, then we force it to 'categorical', because that's all it
    ## realy could have been
    if (mean(is.na(out$type)) > 0.5) {
      out <- mutate(out, type='categorical')
    }
  } else {
    ## Make a dummy, 0 row tibble to send back
    out <- sapply(out.cols, function(x) character(), simplify=FALSE)
    out$date_entered <- integer()
    out <- as.data.frame(out, stringsAsFactors=FALSE) %>% as.tbl
  }

  if (!is.null(covariates)) {
    out <- filter(out, variable %in% covariates)
  }

  out %>% set_fds(x)
}

##' Saves custom sample covariates to a FacileDataSet
##'
##' @export
##' @importFrom jsonlite stream_out
##'
##' @param x the \code{FacileDataSet}
##' @param name the variable name of the covariate
##' @param annotation the annotation table of covariate vaues to a
##'   sample-descriptor-like table
##' @param custom_key the custom key (likely userid) for the annotation
##' @param file.prefix Vincent uses this
##' @param sample_filter_criteria optional list of filtering criteria that were
##'   used to drill down into the samples we have the \code{annotatino}
##'   data.frame for
save_custom_sample_covariates <- function(x, name, annotation,
                                          custom_key=Sys.getenv("USER"),
                                          file.prefix="facile",
                                          sample_filter_critera=NULL) {
  stopifnot(is.FacileDataSet(x))
  stopifnot(is.character(name) && length(name) == 1L)
  if (is.null(custom_key)) {
    custom_key <- 'anonymous'
  }
  custom_key <- make.names(custom_key)
  name <- make.names(name)
  if (annotation$variable[1L] != name) {
    annotation <- mutate(annotation, variable=name)
  }
  annotation <- mutate(annotation, date_entered=as.integer(Sys.time()))

  fn <- paste0(file.prefix, '_', custom_key, '_', name, '_', Sys.Date(),'.json')
  fn <- file.path(x$anno.dir, fn)
  ## TODO: figure out how to encode the sample_filter_criteria into the JSON file
  stream_out(x=annotation, con=file(fn))
  invisible(set_fds(annotation, x))
}

##' Appends covariate columns to a query result
##'
##' Note that this function will force the collection of \code{x}
##'
##' @export
##' @param x a facile sample descriptor
##' @param covariates character vector of covariate names. If \code{NULL}
##'   (default), returns all covariates, if is character and length() == 0, then
##'   this is a no-op (x is returned)
##' @param custom_key The key to use to fetch more custom annotations over
##'   the given samples
##' @param .fds A \code{FacileDataSet} object
##' @return The facile \code{x} object, annotated with the specified covariates.
with_sample_covariates <- function(x, covariates=NULL,
                                   custom_key=Sys.getenv("USER"), .fds=fds(x)) {
  assert_sample_subset(x)
  stopifnot(is.FacileDataSet(.fds))
  stopifnot(is.character(covariates) || is.null(covariates))

  if (is.character(covariates) && length(covariates) == 0L) {
    return(x)
  }

  samples <- x %>%
    select(dataset, sample_id) %>%
    collect(n=Inf) %>% ## can't call distinct on SQLite backend :-(
    distinct(.keep_all=TRUE)

  covs <- fetch_sample_covariates(.fds, samples, covariates,
                                  custom_key=custom_key) %>%
    spread_covariates(.fds)

  # if (is.data.table(samples)) {
  #   setDT(covs)
  # }

  collect(x, n=Inf) %>%
    left_join(covs, by=c('dataset', 'sample_id')) %>%
    set_fds(.fds)
}

##' Spreads the covariates returned from database into wide data.frame
##'
##' Samples that did not have a value for a specific covariate are assigned to
##' have NA.
##'
##' @export
##' @param x output from \code{fetch_sample_covariates}
##' @param .fds A \code{FacileDataSet} object
##' @return a wide \code{tbl_df}-like object
spread_covariates <- function(x, .fds=fds(x)) {
  stopifnot(is.FacileDataSet(.fds))
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

  cov.def <- covariate_definitions(.fds)
  if (!is.null(cov.def)) {
    do.cast <- setdiff(colnames(out), c('dataset', 'sample_id'))
    ## Don't decode categorical variables of type 'user_annotation'
    user.anno <- filter(x, class == 'user_annotation' & type == 'categorical')
    if (nrow(user.anno)) {
      do.cast <- setdiff(do.cast, unique(user.anno$variable))
    }

    for (cname in do.cast) {
      casted <- cast_covariate(cname, out[[cname]], cov.def)
      if (is.data.frame(casted)) {
        ## casting a survival covariate will return a two column thing with time
        ## and censoring information, so we need to account for that.
        out[[cname]] <- NULL
        out <- bind_cols(out, casted)
      } else {
        out[[cname]] <- casted
      }
    }
  }

  set_fds(out, .fds)
}

##' Casts the character values of the covariates to their defined types.
##'
##' For most things, a single value will be returned from each cast, but in the
##' case of "time_to_event" data, the value is expended to a two column
##' data.frame with a \code{tte_<covariate>} column for time to event, and an
##' \code{event_<covariate>} column to indicate event (1) or right censored (2).
##'
##' @export
##' @param covariate the name of the covariate
##' @param values the covariate values (which is a \code{character}) as it is
##'   pulled from the database.
##' @param cov.def the un-yamled covariate definitions, if missing we rely on
##'   pulling this out from the \code{FacileDataSet} object \code{.fds}
##' @param .fds If \code{missing(cov.def)}, this is the \code{FacileDataSet} to
##'   get the covariate definitions from.
##' @return values cast to appropriate type if a valid definition was found for
##'   \code{covariate}, otherwise values is returned "as is". Most of the time
##'   this is a single vector, but others it can be a data.frame (for
##'   \code{right_censored} data, for instance)
cast_covariate <- function(covariate, values, cov.def, .fds) {
  if (missing(cov.def)) {
    stopifnot(is.FacileDataSet(.fds))
    cov.def <- covariate_definitions(.fds)
  }
  stopifnot(is(cov.def, 'CovariateDefinitions'))
  stopifnot(is.character(values))
  stopifnot(is.character(covariate) && length(covariate) == 1L)

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
