##' Fetch data from single assay of choice
##'
##' @export
##' @importFrom rhdf5 h5read
##' @importFrom multiGSEA eigenWeightedMean
##' @inheritParams assay_feature_info
##' @param x A \code{FacileDataSet} object.
##' @param features a feature descriptor (data.frame with assay and feature_id
##'   columms)
##' @param samples a sample descriptor to specify which samples to return data
##'   from.
##' @param normalized return normalize or raw data values, defaults to
##'   \code{raw}
##' @param as.matrix by default, the data is returned in a long-form tbl-like
##'   result. If set to \code{TRUE}, the data is returned as a matrix.
##' @param ... parameters to pass to normalization methods
##' @param subset.threshold sometimes fetching all the genes is faster than
##'   trying to subset. We have to figure out why that is, but I've previously
##'   tested random features of different lengths, and around 700 features was
##'   the elbow.
##' @param aggregate.by do you want individual level results or geneset
##'   scores? Use 'ewm' for eigenWeightedMean, and that's all.
##' @return A lazy \code{\link[dplyr]{tbl}} object with the expression
##'   data to be \code{\link[dplyr]{collect}}ed when \code{db} is provided,
##'   othwerise a \code{tbl_df} of the results.
fetch_assay_data <- function(x, features, samples=NULL, assay_name=NULL,
                             normalized=FALSE, as.matrix=FALSE, ...,
                             subset.threshold=700, aggregate.by=NULL) {
  stopifnot(is.FacileDataSet(x))
  assert_flag(as.matrix)
  assert_flag(normalized)
  assert_number(subset.threshold)

  if (!is.null(assay_name) || is.character(features)) {
    assert_string(assay_name)
    assert_choice(assay_name, assay_names(x))
  }

  if (missing(features)) {
    assert_string(assay_name)
    features <- assay_feature_info(x, assay_name) %>% collect(n=Inf)
  } else {
    if (is.character(features)) {
      features <- tibble(feature_id=features, assay=assay_name)
    }
    stopifnot(is(features, 'tbl') || is(features, 'data.frame'))
    if (is.null(features$assay)) {
      features$assay <- assay_name
    }
    assert_assay_feature_descriptor(features)
  }

  ## Adding check for a 0 row data.frame, because there is some chain of
  ## reactivity that fires in FacileExplorer upon FDS switching that triggers
  ## this when the previous dataset has sample filters entered. This acts
  ## as a defense to that, and also works to handle this strange case from the
  ## backend side, too -- perhaps a user will stumble on this in their analyses?
  if (is.null(samples) || (is.data.frame(samples) && nrow(samples) == 0L)) {
    samples <- samples(x)
  }
  assert_sample_subset(samples)

  assays <- unique(features$assay)
  n.assays <- length(assays)
  if (n.assays > 1L && as.matrix) {
    stop("Fetching from multiple assays requires return in melted form")
  }

  if (!is.null(aggregate.by)) {
    assert_string(aggregate.by)
    assert_choice(aggregate.by, 'ewm')
    stopifnot(n.assays == 1L)
    if (!normalized) {
      warning("You probably don't want to aggregate.by on unnormalized data",
              immediate.=TRUE)
    }
  }

  out <- lapply(assays, function(a) {
    f <- filter(features, assay == a)
    .fetch_assay_data(x, a, f$feature_id, samples, normalized, as.matrix,
                      subset.threshold, aggregate.by)
  })

  if (length(out) == 1L) {
    out <- out[[1L]]
  } else if (!as.matrix) {
    ## We stop if we are asking for a matrix across multiple assays, but maybe
    ## we don't have to ... (in the future, I mean)
    out <- bind_rows(out)
  }

  out
}

.fetch_assay_data <- function(x, assay_name, feature_ids, samples,
                              normalized=FALSE, as.matrix=FALSE,
                              subset.threshold=700, aggregate.by=NULL, ...) {
  stopifnot(is.FacileDataSet(x))
  assert_string(assay_name)
  assert_character(feature_ids, min.len=1L)
  samples <- assert_sample_subset(samples)
  assert_flag(normalized)
  assert_flag(as.matrix)
  assert_number(subset.threshold)
  if (!is.null(aggregate.by)) {
    assert_string(aggregate.by)
    assert_choice(aggregate.by, 'ewm')
  }

  finfo <- assay_feature_info(x, assay_name, feature_ids=feature_ids) %>%
    collect(n=Inf) %>%
    arrange(hdf5_index)
  atype <- finfo$assay_type[1L]
  ftype <- finfo$feature_type[1L]
  sinfo <- assay_sample_info(x, assay_name, samples) %>%
    mutate(samid=paste(dataset, sample_id, sep="__"))
  bad.samples <- is.na(sinfo$hdf5_index)
  if (any(bad.samples)) {
    warning(sum(bad.samples), " samples not found in `", assay_name, "`assay.",
            immediate.=TRUE)
    sinfo <- sinfo[!bad.samples,,drop=FALSE]
  }

  ## DEBUG: Tune chunk size?
  ## As the number of genes you are fetching increases, only subsetting
  ## out a few of them intead of first loading the whole matrix gives you little
  ## value, for instance, over all BRCA tumors (994) these are some timings for
  ## different numbers of genes:
  ##
  ##     Ngenes                   Time (seconds)
  ##     10                       0.5s
  ##     100                      0.8s
  ##     250                      1.2s
  ##     500                      2.6s
  ##     750                      6s seconds
  ##     3000                     112 seconds!
  ##     unpspecified (all 26.5k) 7 seconds!
  ##
  ## I'm using `ridx` as a hack downstream to run around the issue of slowing
  ## down the data by trying to subset many rows via hdf5 instead of loading
  ## then subsetting after (this is so weird)
  ##
  ## TODO: setup unit tests to ensure that ridx subsetting and remapping back
  ## to original genes works
  fetch.some <- nrow(finfo) < subset.threshold
  ridx <- if (fetch.some) finfo$hdf5_index else NULL

  dat <- sinfo %>%
    group_by(dataset) %>%
    do(res={
      ds <- .$dataset[1L]
      hd5.name <- paste('assay', assay_name, ds, sep='/')
      vals <- h5read(hdf5fn(x), hd5.name, list(ridx, .$hdf5_index))
      if (is.null(ridx)) {
        vals <- vals[finfo$hdf5_index,,drop=FALSE]
      }
      dimnames(vals) <- list(finfo$feature_id, .$samid)
      if (normalized) {
        vals <- normalize.assay.matrix(vals, finfo, ., ...)
      }
      vals
    }) %>%
    ungroup

  ## NOTE: We can avoid the monster matrix creation if we only want !as.matrix
  ## returns, but this makes the code easier to reason. We can come back to this
  ## to optimize for speed later. The problem is introduced when the
  ## aggregate.by parameter was introduced
  vals <- do.call(cbind, dat$res)
  # browser()
  if (nrow(vals) == 1L) {
    if (!is.null(aggregate.by)) {
      warning("No assay feature aggregation performed over single feature",
              immediate.=TRUE)
    }
    aggregate.by <- NULL
  }

  if (!is.null(aggregate.by)) {
    scores <- eigenWeightedMean(vals)$score
    vals <- matrix(scores, nrow=1, dimnames=list('score', names(scores)))
  }

  if (!as.matrix) {
    vals <- .melt.assay.matrix(vals, assay_name, atype, ftype, finfo)
    if (!is.null(aggregate.by)) {
      vals[, feature_type := 'aggregated']
      vals[, feature_id := 'aggregated']
      vals[, feature_name := 'aggregated']
    }
    vals <- as.tbl(setDF(vals))
  }

  class(vals) <- c('FacileExpression', class(vals))
  set_fds(vals, x)
}

.melt.assay.matrix <- function(vals, assay_name, atype, ftype, finfo) {
  vals <- as.data.table(vals, keep.rownames=TRUE)
  vals <- melt.data.table(vals, id.vars='rn', variable.factor=FALSE,
                          variable.name='sample_id')
  setnames(vals, 1L, 'feature_id')

  vals[, dataset := sub('__.*$', '', sample_id)]
  vals[, sample_id := sub('^.*__', '', sample_id)]
  vals[, assay := assay_name]
  vals[, assay_type := atype]
  vals[, feature_type := ftype]
  xref <- match(vals$feature_id, finfo$feature_id)
  vals[, feature_name := finfo$name[xref]]
  corder <- c('dataset', 'sample_id', 'assay', 'assay_type',
              'feature_type', 'feature_id', 'feature_name', 'value')
  setcolorder(vals, corder)
  vals
}

##' Helper function to get sample assay data from single or aggregate features
##' @export
fetch_assay_score <- function(x, features, samples=NULL, assay_name=NULL,
                              as.matrix=FALSE, ..., subset.threshold=700) {
  if (is.null(assay_name)) {
    assay_name <- features$assay
  }
  stopifnot(is.character(assay_name), length(unique(asssay_name)) == 1L)
  dat <- fetch_assay_data(x, features, samples=samples, assay_name=NULL,
                          as.matrix=TRUE, normalized=TRUE,
                          subset.threshold=subset.threshold)
  if (nrow(dat) > 1) {
    dat <- matrix(eigenWeightedMean(dat)$score, nrow=1)
  }

}
##' @export
assay_types <- function(x) {
  stopifnot(is.FacileDataSet(x))
  assay_info_tbl(x) %>% collect(n=Inf) %$% assay_type
}

##' @export
assay_names <- function(x) {
  stopifnot(is.FacileDataSet(x))
  assay_info_tbl(x) %>% collect %$% assay
}

##' @export
has_assay <- function(x, assay_name) {
  stopifnot(is.FacileDataSet(x))
  assert_character(assay_name)
  assay_name %in% assay_names(x)
}

##' Utility functions to get row and column indices of rnaseq hdf5 files.
##'
##' This is called to get things like hdf5_index and scaling factors for
##' the samples in a given assay.
##'
##' @export
##' @param x \code{FacileDataSet}
##' @param assay_name the name of the assay
##' @param samples a sample descriptor
##' @return an updated version of \code{samples} decorated with hd5_index,
##'   scaling factors, etc. Note that rows in \code{samples} that do not appear
##'   in \code{assay_name} will be returnd here with NA values for hd5_index and
##'   such.
assay_sample_info <- function(x, assay_name, samples=NULL) {
  stopifnot(is.FacileDataSet(x))
  if (!is.null(samples)) {
    samples <- assert_sample_subset(samples) %>%
      distinct(dataset, sample_id) %>%
      collect(n=Inf)
  }
  feature.type <- assay_feature_type(x, assay_name) ## validate assay_name

  asi <- assay_sample_info_tbl(x) %>%
    filter(assay == assay_name) %>%
    collect(n=Inf)

  if (is.null(samples)) {
    samples <- asi
  } else {
    samples <- left_join(samples, asi, by=c('dataset', 'sample_id'))
  }

  samples
}


##' Returns the feature_type for a given assay
##'
##' The elements of the rows for a given assay all correspond to a particular
##' feature space (ie. feature_type='entrez')
##'
##' @export
##' @param x \code{FacileDataSet}
##' @param assay_name the name of the assay
assay_feature_type <- function(x, assay_name) {
  stopifnot(is.FacileDataSet(x))
  assert_string(assay_name)
  assert_choice(assay_name, assay_names(x))
  assay_info_tbl(x) %>%
    filter(assay == assay_name) %>%
    collect %$%
    feature_type
}

##' Materializes a table with all feature information for a given assay.
##'
##' DEBUG: This logic is unnecessarily complex because I make sure to collect
##' all tables from the database as opposed to copying external tables in and
##' doing an inner_join in the database. I'm doing this becuase we are getting
##' name collections on some of the temporary tables. We get erros like:
##'     Warning: Error in : Table pkdtpohpsu already exists.
##'
##' This fetches the hdf5_index for the assays as well
##' @export
##' @inheritParams assay_feature_type
##' @param feature_ids a character vector of feature_ids
##' @return a \code{tbl_sqlite} result with the feature information for the
##'   features in a specified assay. This
assay_feature_info <- function(x, assay_name, feature_ids=NULL) {
  ## NOTE: This is currently limited to a single assay
  ftype <- assay_feature_type(x, assay_name)
  if (!is.null(feature_ids)) {
    assert_character(feature_ids)
    if (length(feature_ids) == 0) {
      feature_ids <- NULL
    } else {
      feature_ids <- tibble(feature_id=feature_ids)
    }
  }

  afinfo <- assay_feature_info_tbl(x) %>%
    filter(assay == assay_name) %>%
    collect(n=Inf)

  if (!is.null(feature_ids)) {
    afinfo <- inner_join(afinfo, feature_ids, by=c('feature_id'))
  }

  assay.info <- assay_info_tbl(x) %>%
    select(assay, assay_type, feature_type) %>%
    filter(assay == assay_name) %>%
    collect(n=Inf)

  out <- afinfo %>% inner_join(assay.info, by='assay')
  ftype <- out$feature_type[1L]
  finfo <- feature_info_tbl(x)
  if (length(ftype == 1L)) {
    finfo <- filter(finfo, feature_type == ftype)
  } else {
    finfo <- filter(finfo, feature_type %in% ftype)
  }
  finfo <- collect(finfo, n=Inf)

  out %>%
    inner_join(finfo, by=c('feature_type', 'feature_id')) %>%
    set_fds(x)
}

##' Identify the number of each assay run across specific samples
##'
##' @export
##' @param x FacileDataSet
##' @param samples sample descriptor
##' @param with_count return the number of samples in \code{samples} that are
##'   assayed over each assay as a column in \code{return}
##' @return rows from assay_info_tbl that correspond to the assays defined
##'   over the given samples. If no assays are defined over these samples,
##'   you're going to get an empty tibble.
assay_info_over_samples <- function(x, samples, with_count=TRUE) {
  stopifnot(is.FacileDataSet(x))
  assert_sample_subset(samples)
  assays <- assay_sample_info_tbl(x) %>%
    collect %>%
    inner_join(samples, by=c('dataset', 'sample_id'))
  if (with_count) {
    assays <- group_by(assays, assay) %>%
      summarize(n=n()) %>%
      ungroup
  } else {
    assays <- distinct(assays, assay) %>% mutate(n=-1L)
  }

  inner_join(assays, assay_info_tbl(x), by='assay', copy=TRUE) %>%
    arrange(assay)
}



## helper functino to fetch_assay_data
normalize.assay.matrix <- function(vals, feature.info, sample.info,
                                   log=TRUE, prior.count=5, ...) {
  stopifnot(
    nrow(vals) == nrow(feature.info),
    all(rownames(vals) == feature.info$feature_id),
    ncol(vals) == nrow(sample.info),
    all(colnames(vals) == sample.info$samid),
    is.character(feature.info$assay_type),
    length(unique(feature.info$assay_type)) == 1L,
    is.numeric(sample.info$libsize), is.numeric(sample.info$normfactor))
  atype <- feature.info$assay_type[1L]
  libsize <- sample.info$libsize * sample.info$normfactor
  if (atype == 'rnaseq') {
    out <- edgeR::cpm(vals, libsize, log=log, prior.count=prior.count)
  } else {
    warning("No normalization procedure for ", atype, " assay", immediate.=TRUE)
    out <- vals
  }
  out
}

create_assay_feature_descriptor <- function(x, features, assay_name=NULL) {
  ## TODO: Refactor the code inside `fetch_assay_data` to use this.
  stopifnot(is.FacileDataSet(x))

  if (is.character(features) || missing(features)) {
    assert_string(assay_name)
    assert_choice(assay_name, assay_names(x))
  }

  if (missing(features)) {
    features <- assay_feature_info(x, assay_name) %>% collect(n=Inf)
  } else if (is.character(features)) {
    features <- tibble(feature_id=features, assay=assay_name)
  } else if (is.data.frame(features) && is.character(assay_name)) {
    features[['assay']] <- assay_name
  }
  features <- collect(features, n=Inf) ## not necessary?
  assert_assay_feature_descriptor(features, x)
  features
}

##' Append expression values to sample-descriptor
##'
##' @export
##' @param x a samples descriptor
##' @param feature_ids character vector of feature_ids
##' @param with_symbols Do you want gene symbols returned, too?
##' @param .fds A \code{FacileDataSet} object
##' @return a tbl-like result
with_assay_data <- function(samples, features, assay_name=NULL,
                            normalized=FALSE, aggregate.by=NULL, spread=NULL,
                            ..., .fds=fds(samples)) {
  stopifnot(is.FacileDataSet(.fds))
  assert_sample_subset(samples)
  assert_flag(normalized)
  if (!is.null(spread)) {
    spread <- assert_choice(spread, c('id', 'name'))
    spread <- if (spread == 'id') 'feature_id' else 'feature_name'
  }
  features <- create_assay_feature_descriptor(.fds, features, assay_name=assay_name)
  assays <- unique(features$assay)
  if (length(assays) > 1L && !is.null(spread)) {
    stop("Can only spread assay_data if asking for one assay_name")
  }
  adata <- fetch_assay_data(.fds, features, samples, normalized=normalized,
                            aggregate.by=aggregate.by)
  if (!is.null(spread)) {
    adata <- select_(adata, .dots=c('dataset', 'sample_id', spread, 'value'))
    adata <- tidyr::spread_(adata, spread, 'value') %>% set_fds(.fds)
  }

  adata %>% join_samples(samples)
}

##' Takes a result from fetch_expression and spreads out genes acorss columns
##'
##' This is a convenience function, and will try to guess what you mean if you
##' don't explicitly specify which columns to spread and what to call them.
##' With that mind set, if we find a cpm and symbol column, we will use them
##' because those are the thing you will likely want to use for exploratory
##' data analysis if they're in the incoming dataset. If those columns aren't
##' found, then we'll pick the feature_id and count column.
##'
##' @export
##' @param x facile expression result from \code{fetch_expression}
##' @param key the column from the long-form \code{fetch_expression} table
##'   to put in the columns of the outgoing data.frame that the values are
##'   "spread into"
##' @param value the value column to spread into the \code{key} columns
##' @param .fds the \code{FacileDataSet}
##' @return a more stout \code{x} with the expression values spread across
##'   columns.
spread_assay_data <- function(x, assay_name, key=c('symbol', 'feature_id'),
                              .fds=fds(x)) {
  stop("Put spread argument in with_assay_data (is this right?)")
  force(.fds)
  if (missing(key)) {
    key <- if ('symbol' %in% colnames(x)) 'symbol' else 'feature_id'
  }
  if (missing(value)) {
    value <- if ('cpm' %in% colnames(x)) 'cpm' else 'count'
  }
  key <- match.arg(key, c('symbol', 'feature_id'))
  value <- match.arg(value, c('count', 'cpm'))
  assert_expression_result(x)
  assert_columns(x, c(key, value))
  x <- collect(x, n=Inf)

  if (key == 'symbol') {
    f2s <- distinct(x, feature_id, symbol)
    if (any(is.na(f2s$symbol))) stop("NAs found in symbol column, spread")
    if (any(duplicated(f2s$symbol))) stop("Duplicate symbols found")
    x <- select(x, -feature_id)
  } else {
    if ('symbol' %in% colnames(x)) {
      x <- select(x, -symbol)
    }
    x <- mutate(x, feature_id=paste0('feature_id_', feature_id))
  }

  if (value == 'cpm') {
    x <- select(x, -count)
  } else if ('cpm' %in% colnames(x)) {
    x <- select(x, -cpm)
  }

  out <- spread_(x, key, value) %>% set_fds(.fds)
  if (nrow(out) >= nrow(x)) {
    ## You might be tempted to test the width of the outgoing object, too, but
    ## if you only had to features in this object, then it wouldn't have changed
    stop("The spread_operation did not make your incoming object more stout, ",
         "You need to debug this")
  }
  out
}
