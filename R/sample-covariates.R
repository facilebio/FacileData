#' Fetch rows from sample_covariate table for specified samples and covariates
#'
#' @export
#' @param db a \code{FacileDataSet} connection
#' @param samples a samples descriptor \code{tbl_*}
#' @param covariates character vector of covariate names
#' @param custom_key The key to use to fetch more custom annotations over
#'   the given samples
#' @return rows from the \code{sample_covariate} table
#' @family API
fetch_sample_covariates.FacileDataSet <- function(x, samples=NULL, covariates=NULL,
                                    custom_key=Sys.getenv("USER")) {
  stopifnot(is.FacileDataSet(x))
  ## db temp table thing shouldn't be an issue here
  # dat <- sample_covariate_tbl(x) %>% collect(n=Inf) ## #dboptimize# remove to exercise db harder
  dat <- sample_covariate_tbl(x)
  if (is.character(covariates)) {
    if (length(covariates) == 1L) {
      dat <- filter(dat, variable == covariates)
    } else if (length(covariates) > 1L) {
      dat <- filter(dat, variable %in% covariates)
    }
  }
  dat <- collect(dat, n=Inf)
  dat <- set_fds(dat, x) ## explicitly added here to do `collect` above

  ## If the samples descriptor is defined over the sample_covariate table,
  ## this thing explodes (inner joining within itself, I guess). We defensively
  ## copy the sample descriptor, but in future maybe better to test if the
  ## dat and samples sqlite tables are pointing to the same thing
  if (!is.null(samples)) {
    samples <- assert_sample_subset(samples) %>%
      distinct(dataset, sample_id) %>%
      collect(n=Inf)
  }
  # out <- filter_samples(dat, samples)
  out <- join_samples(dat, samples, semi=TRUE)

  if (!is.null(custom_key)) {
    custom <- fetch_custom_sample_covariates(x, samples,
                                             covariates=covariates,
                                             custom_key)
    out <- bind_rows(collect(out, n=Inf), custom)
  }

  set_fds(out, x)
}

#' Fetches custom (user) annotations for a given user prefix
#'
#' @export
#' @importFrom jsonlite stream_in
#' @param fds The \code{FacileDataSet}
#' @param samples the facile sample descriptor
#' @param custom_key The key to use for the custom annotation
#' @return covariate tbl
#' @family API
fetch_custom_sample_covariates.FacileDataSet <- function(x, samples=NULL, covariates=NULL,
                                           custom_key=Sys.getenv("USER"),
                                           file.prefix="facile") {
  stopifnot(is.FacileDataSet(x))
  out.cols <- colnames(sample_covariate_tbl(x))

  fpat <- paste0('^', file.prefix, '_', custom_key, "_.*json")
  annot.files <- list.files(path=x$anno.dir, pattern=fpat, full.names=TRUE)

  if (length(annot.files)) {
    annos <- lapply(annot.files, function(fn) stream_in(file(fn), verbose=FALSE))
    out <- bind_rows(annos) %>%
      select_(.dots=out.cols) %>%
      set_fds(x) %>%
      # filter_samples(samples)
      join_samples(samples, semi=TRUE)
    ## We weren't saving the type == 'categorical' column earlier. So if this
    ## column is.na, then we force it to 'categorical', because that's all it
    ## realy could have been
    if (nrow(out) && mean(is.na(out$type)) > 0.5) {
      out <- mutate(out, class='categorical', type='user_annotation')
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

#' Saves custom sample covariates to a FacileDataSet
#'
#' @export
#' @importFrom jsonlite stream_out
#'
#' @param x the \code{FacileDataSet}
#' @param annotation the annotation table of covariate values to a
#'   sample-descriptor-like table
#' @param name the variable name of the covariate
#' @param custom_key the custom key (likely userid) for the annotation
#' @param file.prefix Vincent uses this
#' @param sample_filter_criteria optional list of filtering criteria that were
#'   used to drill down into the samples we have the \code{annotatino}
#'   data.frame for
save_custom_sample_covariates <- function(x, annotation, name=NULL,
                                          class='categorical',
                                          custom_key=Sys.getenv("USER"),
                                          file.prefix="facile",
                                          sample_filter_critera=NULL) {
  #' TODO: Figure out how to encode sample_filter_criteria into serialized
  #' (JSON) annotation file
  stopifnot(is.FacileDataSet(x))
  annotation <- collect(annotation, n=Inf)
  assert_columns(annotation, c('dataset', 'sample_id', 'value'))
  if (is.null(name)) name <- annotation$name
  if (!test_string(name)) stop("No name given/inferred for custom annotation")

  if (is.null(custom_key)) custom_key <- 'anonymous'
  custom_key <- assert_string(custom_key) %>% make.names

  annotation[['variable']] <- make.names(name)
  annotation <- annotation[, c('dataset', 'sample_id', 'variable', 'value')]
  annotation[['class']] <- class
  annotation[['type']] <- 'user_annotation'
  annotation[['date_entered']] <- as.integer(Sys.time())

  fn <- paste0(file.prefix, '_', custom_key, '_', name, '_', Sys.Date(),'.json')
  fn <- file.path(x$anno.dir, fn)
  stream_out(x=annotation, con=file(fn))
  invisible(set_fds(annotation, x))
}

#' Appends covariate columns to a query result
#'
#' Note that this function will force the collection of \code{x}
#'
#' @export
#' @importFrom stats complete.cases
#' @param x a facile sample descriptor
#' @param covariates character vector of covariate names. If \code{NULL}
#'   (default), returns all covariates, if is character and length() == 0, then
#'   this is a no-op (x is returned)
#' @param na.rm if \code{TRUE}, filters outgoing result such that only rows
#'   with nonNA values for the \code{covariates} specified here will be
#'   returned. Default: \code{FALSE}. Note that this will not check columns
#'   not specified in \code{covariates} for NA-ness.
#' @param custom_key The key to use to fetch more custom annotations over
#'   the given samples
#' @param .fds A \code{FacileDataSet} object
#' @return The facile \code{x} object, annotated with the specified covariates.
with_sample_covariates <- function(x, covariates=NULL, na.rm=FALSE,
                                   custom_key=Sys.getenv("USER"), .fds=fds(x)) {
  stopifnot(is.FacileDataSet(.fds))
  x <- assert_sample_subset(x) %>% collect(n=Inf)
  stopifnot(is.character(covariates) || is.null(covariates))
  if (is.character(covariates) && length(covariates) == 0L) {
    return(x)
  }

  samples <- x %>%
    select(dataset, sample_id) %>%
    distinct(.keep_all=TRUE)

  covs <- fetch_sample_covariates(.fds, samples, covariates,
                                  custom_key=custom_key)
  covs <- spread_covariates(covs, .fds)

  out <- left_join(x, covs, by=c("dataset", "sample_id"), suffix = c("", ".x"))

  if (na.rm && length(covariates)) {
    keep <- complete.cases(select_(out, .dots=covariates))
    out <- out[keep,,drop=FALSE]
  }

  set_fds(out, .fds)
}

#' Spreads the covariates returned from database into wide data.frame
#'
#' Samples that did not have a value for a specific covariate are assigned to
#' have NA.
#'
#' @export
#' @param x output from \code{fetch_sample_covariates}
#' @param .fds A \code{FacileDataSet} object
#' @return a wide \code{tbl_df}-like object
spread_covariates <- function(x, .fds=fds(x)) {
  stopifnot(is.FacileDataSet(.fds))
  x <- assert_sample_covariates(x) %>%
    collect(n=Inf)

  ## Ensures we get a row for every sample in x, even if it is missing a value
  ## for the covariate
  dummy <- select(x, dataset, sample_id) %>%
    collect(n=Inf) %>%
    distinct(.keep_all=TRUE) %>%
    mutate(variable='.dummy.', value=NA)

  out <- bind_rows(x, dummy) %>%
    dcast(dataset + sample_id ~ variable, value.var='value') %>%
    mutate(.dummy.=NULL) %>%
    ## I don't think we should set rownames here. Delete next command and test
    # set_rownames(., paste(.$dataset, .$sample_id, sep='__'))
    as.tbl

  cov.def <- covariate_definitions(.fds)
  if (!is.null(cov.def)) {
    do.cast <- setdiff(colnames(out), c('dataset', 'sample_id'))
    ## Don't decode categorical variables of type 'user_annotation'
    user.anno <- filter(x, type == 'user_annotation' & class == 'categorical')
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
