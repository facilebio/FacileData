#' @export
#' @noRd
fetch_assay_data.facile_frame <- function(x, features = NULL, samples = NULL,
                                          assay_name = NULL,
                                          normalized = FALSE,
                                          batch = NULL, main = NULL,
                                          as.matrix = FALSE,
                                          drop_samples = TRUE,
                                          ...,
                                          subset.threshold = 700,
                                          aggregate = FALSE,
                                          aggregate.by= "ewm",
                                          verbose = FALSE) {
  fds. <- assert_facile_data_store(fds(x))
  if (!is.null(samples)) {
    warning("`samples` ignored when fetching assay_data from a facile_frame",
            immediate. = TRUE)
  }
  if (is.null(assay_name)) {
    assay_name <- default_assay(fds.)
  }

  fetch_assay_data(fds., features = features, samples = x,
                   assay_name = assay_name, normalized = normalized,
                   batch = batch, main = main,
                   as.matrix = as.matrix, drop_samples = drop_samples, ...,
                   subset.threshold = subset.threshold,
                   aggregate = aggregate, aggregate.by = aggregate.by,
                   verbose = verbose)
}

#' @noRd
#' @export
assay_names.facile_frame <- function(x, default_first = TRUE, ...) {
  fds. <- fds(x)
  all_names <- assay_names(fds., default_first = default_first, ...)
  anames <- lapply(all_names, \(aname) {
    asi <- assay_sample_info(x, aname, drop_samples = TRUE)
    if (nrow(asi)) aname else character()
  })
  unlist(anames, use.names = FALSE)
}

#' @export
#' @noRd
assay_info.FacileDataSet <- function(x, assay_name = NULL, ...) {
  stopifnot(is.FacileDataSet(x))
  ainfo <- assay_info_tbl(x) |> collect(n = Inf)
  if (!is.null(assay_name)) {
    assert_string(assay_name)
    assert_choice(assay_name, ainfo$assay)
    ainfo <- filter(ainfo, assay == assay_name)
  }
  as_facile_frame(ainfo, x)
}

#' @noRd
#' @export
assay_sample_info.facile_frame <- function(x, assay_name, drop_samples = TRUE,
                                           ...) {
  assay_sample_info(fds(x), assay_name, drop_samples = drop_samples,
                    samples = collect(x, n = Inf), ...)
}

#' Returns the feature_type for a given assay
#'
#' The elements of the rows for a given assay all correspond to a particular
#' feature space (ie. feature_type='entrez')
#'
#' @export
#' @param x A `FacileDataStore`
#' @param assay_name the name of the assay
assay_feature_type <- function(x, assay_name) {
  assert_class(x, "FacileDataStore")
  assert_choice(assay_name, assay_names(x))
  assay_info(x) |>
    filter(.data$assay == .env$assay_name) |>
    collect(n = Inf) |>
    pull(feature_type)
}


#' @rdname feature_name_map
#' @export
#'
#' @param x \code{FacileDataSet}
#' @param assay_name the name of assay to get the feature map for.
assay_feature_name_map <- function(x, assay_name) {
  ftype <- assay_feature_type(x, assay_name)
  feature_name_map(x, ftype)
}


#' Identify the number of each assay run across specific samples.
#'
#' The default assay is listed first, the rest of the order is undetermined.
#'
#' ```
#' A tibble: 2 × 3
#' assay    ndatasets nsamples
#' <chr>        <int>    <int>
#' scrnaseq         3       26
#' snrnaseq         3       20
#' ```
#'
#' @export
#' @param x FacileDataSet
#' @param samples sample descriptor
#' @param with_count return the number of samples in \code{samples} that are
#'   assayed over each assay as a column in \code{return}
#' @return rows from assay_info_tbl that correspond to the assays defined
#'   over the given samples. If no assays are defined over these samples,
#'   you're going to get an empty tibble.
assay_summary <- function(x, ...) {
  UseMethod("assay_summary", x)
}

#' @noRd
#' @export
assay_summary.FacileDataStore <- function(x, samples = NULL, ...) {
  if (!is.null(samples)) {
    assert_sample_subset(samples, x)
    samples <- distinct(samples, dataset, sample_id)
  } else {
    samples <- collect(samples(x), n = Inf)
  }
  assay_summary(samples, ...)
}  

#' @noRd
#' @export
assay_summary.facile_frame <- function(x, ...) {
  fds. <- fds(x)
  anames <- assay_names(fds.)
  info <- lapply(anames, function(aname) {
    x |> 
      assay_sample_info(aname, drop_samples = TRUE) |> 
      summarize(ndatasets = n_distinct(dataset), nsamples = n()) |> 
      mutate(assay = aname, .before = 1L)
  })
  
  bind_rows(info) |> 
    filter(.data$nsamples > 0L) |> 
    left_join(
      select(assay_info(fds.), -description, -storage_mode),
      by = "assay")
}

#' Check if character vector of sample ids in `ids` can plausably be the
#' sample_id's in the sample_info table.
#'
#' @noRd
samples_look_concordant <- function(ids, sample_info) {
  assert_multi_class(sample_info, c("data.frame", "tbl"))
  assert_character(ids, len = nrow(sample_info))
  if (is.element("sample_id", colnames(sample_info))) {
    sids <- sample_info[["sample_id"]]
    # let's take the last N entries from ids, where N is the length of the
    # sister entry in sids
    ids.trim <- substr(ids, nchar(ids) - nchar(sids) + 1L, nchar(ids))
    kosher <- isTRUE(all.equal(ids.trim, sids))
  } else {
    sids <- rownames(sample_info)
    kosher <- isTRUE(all.equal(ids, sids))
  }
  kosher
}

.extract_batch_terms <- function(sample.info, batch, batch2, covariates,
                                 design, batch_columns, ...) {
  # terms()
}

#' @export
#' @rdname fetch_assay_data
#'
#' @param samples a samples descriptor
#' @param feature_ids character vector of feature_ids
#' @param with_symbols Do you want gene symbols returned, too?
#' @param .fds A \code{FacileDataSet} object
#' @return a tbl-like result
with_assay_data.facile_frame <- function(x, features, assay_name = NULL,
                                         normalized = TRUE, aggregate = FALSE,
                                         aggregate.by = "ewm", spread = TRUE,
                                         with_assay_name = FALSE, ...,
                                         verbose = FALSE, .fds = fds(x)) {
  x <- collect(x, n = Inf)
  .fds <- assert_facile_data_store(.fds)
  NextMethod(x, .fds = .fds)
}

#' @export
#' @noRd
with_assay_data.tbl <- function(x, features, assay_name = NULL,
                                normalized = TRUE, aggregate = FALSE,
                                aggregate.by = "ewm", spread = TRUE,
                                with_assay_name = FALSE, ..., verbose = FALSE,
                                .fds = NULL) {
  with_assay_data.data.frame(collect(x, n = Inf), features = features,
                             assay_name = assay_name, normalized = normalized,
                             aggregate = aggregate, aggregate.by = aggregate.by,
                             spread = spread, with_assay_name = with_assay_name,
                             ..., verbose = verbose, .fds = .fds)
}

#' @export
#' @method with_assay_data data.frame
#' @noRd
with_assay_data.data.frame <- function(x, features, assay_name = NULL,
                                       normalized = TRUE, aggregate = FALSE,
                                       aggregate.by = "ewm",
                                       spread = TRUE, with_assay_name=FALSE,
                                       ..., verbose = FALSE, .fds = NULL) {
  .fds <- assert_facile_data_store(.fds)
  assert_sample_subset(x)
  # assert_flag(normalized)

  ## Check that parameters are kosher before fetching data
  features <- create_assay_feature_descriptor(.fds, features,
                                              assay_name = assay_name)
  assay_name <- unique(features$assay)
  if (test_flag(spread) && spread) {
    ## infer column based on assay type (rnaseq for now)
    spread <- if (can.spread.assay.by.name(adata, assay_name)) 'name' else 'id'
  }
  if (is.character(spread)) {
    spread <- assert_choice(spread, c('id', 'name'))
    spread <- if (spread == 'id') 'feature_id' else 'feature_name'
  }
  if (length(assay_name) > 1L && !is.null(spread)) {
    stop("Can only spread assay_data if asking for one assay_name")
  }

  # Hit the datastore
  adata <- fetch_assay_data(.fds, features, x, assay_name = assay_name,
                            normalized = normalized, aggregate = aggregate,
                            aggregate.by = aggregate.by, verbose = verbose, ...)
  adata <- filter(adata, !is.na(value))

  if (is.character(spread)) {
    spread.vals <- unique(adata[[spread]])
    if (any(spread.vals %in% colnames(x))) {
      if (!with_assay_name && verbose) {
        warning("appending assay_name to spread columns to avoid collision")
      }
      with_assay_name <- TRUE
    }
    # adata <- select_(adata, .dots=c('dataset', 'sample_id', spread, 'value'))
    adata <- select(adata, dataset, sample_id, !!spread, value)
    adata <- tidyr::spread(adata, spread, 'value')
    spread.idx <- which(colnames(adata) %in% spread.vals)
    if (with_assay_name || spread == 'id') {
      newname <- paste0(assay_name, '_', colnames(adata)[spread.idx])
      colnames(adata)[spread.idx] <- newname
    }
    adata <- set_fds(adata, .fds)
  }

  out <- join_samples(x, adata)
  as_facile_frame(out, .fds)
}

can.spread.assay.by.name <- function(x, assay_name) {
  ## TODO: check if duplicate sample_id;name combos exist, in which case
  ## we spread with id and not name
  TRUE
}

#' Takes a result from fetch_expression and spreads out genes across columns
#'
#' This is a convenience function, and will try to guess what you mean if you
#' don't explicitly specify which columns to spread and what to call them.
#' With that mind set, if we find a cpm and symbol column, we will use them
#' because those are the thing you will likely want to use for exploratory
#' data analysis if they're in the incoming dataset. If those columns aren't
#' found, then we'll pick the feature_id and count column.
#'
#' @export
#' @importFrom stats setNames
#' @param x facile expression result from \code{fetch_expression}
#' @param key the column from the long-form \code{fetch_expression} table
#'   to put in the columns of the outgoing data.frame that the values are
#'   "spread into"
#' @param value the value column to spread into the \code{key} columns
#' @param .fds the \code{FacileDataSet}
#' @return a more stout \code{x} with the expression values spread across
#'   columns.
spread_assay_data <- function(x, assay_name, key=c('name', 'feature_id'),
                              value=c('cpm', 'value', 'count'),
                              .fds=fds(x)) {
  stop("Put spread argument in with_assay_data (is this right?)")
  force(.fds)
  if (missing(key)) {
    key <- if ('name' %in% colnames(x)) 'symbol' else 'feature_id'
  }
  key <- match.arg(key, c('name', 'feature_id'))

  val.opts <- intersect(c('value', 'cpm', 'count'), colnames(x))

  if (missing(value)) {
    xref <- setNames(match(colnames(x), val.opts), colnames(x))
    xref <- xref[!is.na(xref)]

    ## get either of these options in this order
    vals <-
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

  out <- tidyr::spread(x, key, value)
  out <- set_fds(out, .fds)
  if (nrow(out) >= nrow(x)) {
    ## You might be tempted to test the width of the outgoing object, too, but
    ## if you only had to features in this object, then it wouldn't have changed
    stop("The spread_operation did not make your incoming object more stout, ",
         "You need to debug this")
  }
  out
}

#' @noRd
#' @export
with_assay_covariates.facile_frame <- function(x, covariates = NULL,
                                               assay_name = default_assay(x),
                                               ..., .fds = fds(x)) {
  x <- collect(x, n = Inf)
  .fds <- assert_facile_data_store(.fds)

  NextMethod(x, .fds = .fds)
  # NextMethod(x = x, .fds = .fds)
}

#' @noRd
#' @export
#' @method with_assay_covariates data.frame
with_assay_covariates.data.frame <- function(x, covariates = NULL,
                                             assay_name = NULL, ...,
                                             .fds = NULL) {
  # doing fds(.fds) here because of the collision between
  # FacileShine::ReactiveFacileDataSet (boxd) objects and Issue #2
  # Currently we are accessing directly the assay_sample_info tbl to get
  # assay_sample covariates, which needs to change.
  #
  # https://github.com/denalitherapeutics/FacileData/issues/2
  ofds <- .fds
  .fds <- assert_facile_data_store(fds(.fds))

  # Until Issue #2 is complete, we can only fetch "libsize" and "normfactor"
  choices <- c("libsize", "normfactor")
  if (is.null(covariates)) covariates <- choices
  assert_character(covariates, null.ok = TRUE)
  assert_subset(covariates, choices, empty.ok = TRUE)

  if (is.null(assay_name)) assay_name <- default_assay(.fds)
  assert_choice(assay_name, assay_names(.fds))


  x <- collect(x, n = Inf)
  assert_sample_subset(x, .fds)

  ss <- assay_sample_info_tbl(.fds) |>
    filter(assay == assay_name) |>
    select(dataset, sample_id, !!covariates) |>
    collect(n = Inf)

  out <- left_join(x, ss, by = c("dataset", "sample_id"))
  as_facile_frame(out, ofds)
}

#' Helper function to get sample assay data from single or aggregate features
#' @export
#' @family API
fetch_assay_score.FacileDataSet <- function(x, features, samples = NULL,
                                            assay_name = NULL,
                                            as.matrix = FALSE, ...,
                                            subset.threshold = 700) {
  .Deprecated("fatch_assay_data(..., aggregate = TRUE)")
  if (is.null(assay_name)) {
    assay_name <- features$assay
  }
  stopifnot(is.character(assay_name), length(unique(asssay_name)) == 1L)
  dat <- fetch_assay_data(x, features, samples = samples, assay_name = NULL,
                          as.matrix = TRUE, normalized = TRUE, ...,
                          subset.threshold = subset.threshold)
  if (nrow(dat) > 1) {
    dat <- matrix(sparrow::eigenWeightedMean(dat)$score, nrow = 1)
  }

}
