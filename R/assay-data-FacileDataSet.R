# These functions are implementations around the facile assay functions that
# are specfici to a FacileDataSet — these should probably be moved to their
# own FacileDataSet package.

#' @export
#' @noRd
fetch_assay_data.FacileDataSet <- function(x, features = NULL, samples = NULL,
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
  assert_flag(as.matrix)
  assert_flag(normalized)
  assert_number(subset.threshold)
  if (!is.null(batch)) {
    normalize <- TRUE
    assert_character(batch, min.len = 1L)
    if (!is.null(main)) {
      assert_character(main)
      if (length(main) == 0L) main <- NULL
    }
  }
  
  features <- create_assay_feature_descriptor(x, features, assay_name = assay_name)
  assay_name <- unique(features$assay)

  assert_string(assay_name)
  assert_choice(assay_name, assay_names(x))
  
  # This was added in 2024-11-22, because I don't know why I thought we'd ask
  # for multiple assay retrieval. These should be different calls to fetch_*
  all.assays <- unique(features$assay)
  if (length(all.assays) != 1L) {
    stop("We should only allow 1 assay per fetch here")
  }
  if (all.assays != assay_name) {
    stop("assay_name is different from features()$assay: ",
         "assay_name: ", assay_name, " <> features()$assay: ", all.assays)
  }

  # I had originally add tests for 0-row datasets returning data from all
  # samples (as if it were NULL) but this was throwing me for a loop in analysis
  # code. If sampels is a 0-row sample-descriptor, the user (like me, ugh) may
  # have inadvertently filtered away all the samples without realizing and query
  # for data and not realize you are getting data from samples you weren't
  # expecting. (what's up with the diary entry?)
  if (is.null(samples)) {
    samples <- samples(x)
    extra_classes <- NULL
  } else {
    ignore.tbl <- class(samples)
    ignore.tbl <- ignore.tbl[grepl("^tbl", ignore.tbl)]
    ignore <- c("data.frame", ignore.tbl)
    extra_classes <- setdiff(class(samples), ignore)
    samples <- as_facile_frame(samples, x) # does validity checks too
  }
  
  samples <- collect(samples, n = Inf)
  # you might think that you want to refactor this and put it before the
  # `collect()` call, but the SQLite back end can't handle
  # `distinct(..., .keep_all = TRUE)`
  samples <- distinct(samples, dataset, sample_id, .keep_all = TRUE)
  
  if (nrow(samples) == 0) {
    warning("Empty sample descriptor provided", immediate. = TRUE)
    return(samples)
  }
  
  n.assays <- length(all.assays)
  if (n.assays > 1L && as.matrix) {
    stop("Fetching from multiple assays requires return in melted form")
  }
  
  if (isTRUE(aggregate)) {
    assert_string(aggregate.by)
    aggregate.by <- assert_choice(tolower(aggregate.by), c('ewm', 'zscore'))
    stopifnot(n.assays == 1L)
    if (!normalized) {
      warning("You probably don't want to aggregate.by on unnormalized data",
              immediate.=TRUE)
    }
  }
  
  out <- sapply(all.assays, function(a) {
    f <- filter(features, .data$assay == .env$a)
    .fetch_assay_data(x, a, f$feature_id, samples, normalized,
                      batch, main, as.matrix, drop_samples,
                      subset.threshold, aggregate, aggregate.by, ...,
                      verbose = verbose)
  }, simplify = FALSE)
  
  ftype <- features$feature_type[1L]
  atype <- features$assay_type[1L]
  
  aggregated.stats <- lapply(out, function(xx) {
    agg <- attr(xx, "aggregated")
    if (is.list(agg)) {
      agg$method <- aggregate.by
      agg$weights <- dplyr::tibble(
        feature_id = names(agg$weights),
        weight = unname(agg$weights)) |> 
        dplyr::arrange(dplyr::desc(weight))
    }
    agg
  })
  aggregated.stats <- aggregated.stats[[1L]]
  dropped.samples <- lapply(out, attr, "samples_dropped") |> dplyr::bind_rows()
  
  if (length(out) == 1L) {
    out <- out[[1L]]
  } else {
    # We alredy stop()ped if we were asked for a matrix across multiple assays,
    # so `out` must be populated with tibbles
    out <- bind_rows(out)
  }
  
  if (!as.matrix) {
    out <- as_facile_frame(out, x, classes = extra_classes,
                           .valid_sample_check = FALSE)
  } else {
    info <- strsplit(colnames(out), "__")
    attr(out, "samples") <- tibble(
      dataset = sapply(info, "[[", 1L),
      sample_id = sapply(info, "[[", 2L))
  }
  
  attr(out, "assay_name") <- assay_name
  attr(out, "feature_type") <- ftype
  attr(out, "assay_type") <- atype
  attr(out, "samples_dropped") <- dropped.samples
  attr(out, "aggregated") <- aggregated.stats
  out
}


#' @noRd
#' @importFrom rhdf5 h5read
#' @importFrom data.table setDF
.fetch_assay_data <- function(x, assay_name, feature_ids, samples,
                              normalized = FALSE, batch = NULL, main = NULL,
                              as.matrix = FALSE, drop_samples = TRUE,
                              subset.threshold = 700, aggregate = FALSE,
                              aggregate.by = "ewm", ...,
                              verbose=FALSE) {
  stopifnot(is.FacileDataSet(x))
  assert_string(assay_name)
  assert_character(feature_ids, min.len=1L)
  samples <- assert_sample_subset(samples)
  assert_flag(normalized)
  assert_flag(as.matrix)
  assert_number(subset.threshold)
  if (isTRUE(aggregate)) {
    assert_string(aggregate.by)
    aggregate.by <- assert_choice(tolower(aggregate.by), c('ewm', 'zscore'))
  }
  
  finfo <- features(x, assay_name, feature_ids = feature_ids) |>
    collect(n = Inf) |>
    arrange(hdf5_index)
  missed <- setdiff(feature_ids, finfo$feature_id)
  if (nrow(finfo) == 0L) {
    stop("No requested features were found")
  }
  if (length(missed) > 0L) {
    msg <- sprintf(
      "%d/%d features not found in assay: ",
      length(missed),
      length(feature_ids))
    warning(msg, immediate. = TRUE)
  }
  
  atype <- finfo$assay_type[1L]
  ftype <- finfo$feature_type[1L]
  sinfo <- assay_sample_info(samples, assay_name, drop_samples = FALSE)
  sinfo <- mutate(sinfo, samid = paste(dataset, sample_id, sep = "__"))
  
  bad.samples <- is.na(sinfo$hdf5_index)
  dropped.samples <- sinfo |> 
    filter(bad.samples) |> 
    select(dataset, sample_id)
  if (nrow(dropped.samples)) {
    if (verbose) {
      warning(sum(bad.samples), " samples not found in `",
              assay_name, "`assay.", immediate.=TRUE)
    }
    sinfo <- sinfo[!bad.samples,,drop=FALSE]
  }
  
  # DEBUG: Tune chunk size?
  # As the number of genes you are fetching increases, only subsetting
  # out a few of them intead of first loading the whole matrix gives you little
  # value, for instance, over all BRCA tumors (994) these are some timings for
  # different numbers of genes:
  #
  #     Ngenes                   Time (seconds)
  #     10                       0.5s
  #     100                      0.8s
  #     250                      1.2s
  #     500                      2.6s
  #     750                      6s seconds
  #     3000                     112 seconds!
  #     unpspecified (all 26.5k) 7 seconds!
  #
  # I'm using `ridx` as a hack downstream to run around the issue of slowing
  # down the data by trying to subset many rows via hdf5 instead of loading
  # then subsetting after (this is so weird)
  #
  # TODO: setup unit tests to ensure that ridx subsetting and remapping back
  # to original genes works
  fetch.some <- nrow(finfo) < subset.threshold
  ridx <- if (fetch.some) finfo$hdf5_index else NULL
  
  dat <- sinfo |>
    group_by(dataset) |>
    do(res = {
      ds <- .$dataset[1L]
      hd5.name <- paste('assay', assay_name, ds, sep='/')
      vals <- h5read(hdf5fn(x), hd5.name, list(ridx, .$hdf5_index))
      if (is.null(ridx)) {
        vals <- vals[finfo$hdf5_index,,drop=FALSE]
      }
      dimnames(vals) <- list(finfo$feature_id, .$samid)
      vals
    }) |>
    ungroup()
  
  # NOTE: We can avoid the monster matrix creation if we only want !as.matrix
  # returns, but this makes the code easier to reason about.
  # We can come back to this to optimize for speed later. The problem is
  # introduced when the aggregate.by parameter was introduced
  vals <- do.call(cbind, dat$res)
  
  if (normalized) {
    xref <- match(colnames(vals), sinfo[["samid"]])
    if (!isTRUE(all.equal(xref, seq(nrow(xref))))) {
      sinfo <- sinfo[xref,,drop = FALSE]
    }
    vals <- normalize_assay_data(vals, finfo, sinfo, batch = batch, main = main,
                                 verbose = verbose, ..., .fds = x)
  }
  
  if (nrow(vals) == 1L) {
    if (isTRUE(aggregate) && verbose) {
      warning("No assay feature aggregation performed over single feature",
              immediate.=TRUE)
    }
    aggregate <- FALSE
  }
  
  aggregated <- NULL
  if (isTRUE(aggregate)) {
    if (aggregate.by == "ewm") {
      aggregated <- sparrow::eigenWeightedMean(vals, ...)
    } else if (aggregate.by == "zscore") {
      aggregated <- sparrow::zScore(vals, ...)
    } else {
      stop("Unknown aggregation method: ", aggregate.by)
    }
    
    # scores <- switch(aggregate.by,
    #                  ewm = sparrow::eigenWeightedMean(vals, ...)$score,
    #                  zscore = sparrow::zScore(vals, ...)$score)
    # vals <- matrix(scores, nrow=1, dimnames=list('score', names(scores)))
    vals <- matrix(
      aggregated$score,
      nrow = 1L,
      dimnames = list("score", names(aggregated$score)))
  }
  
  if (!as.matrix) {
    vals <- .melt.assay.matrix(vals, assay_name, atype, ftype, finfo)
    if (isTRUE(aggregate)) {
      # vals[, feature_type := 'aggregated']
      # vals[, feature_id := 'aggregated']
      # vals[, feature_name := 'aggregated']
      data.table::set(vals, j = "feature_type", value = "aggregated")
      data.table::set(vals, j = "feature_id", value = "aggregated")
      data.table::set(vals, j = "feature_name", value = "aggregated")
    }
    vals <- as_tibble(setDF(vals))
    if (drop_samples) {
      vals <- inner_join(samples, vals, by = c("dataset", "sample_id"))
    } else {
      vals <- left_join(samples, vals, by = c("dataset", "sample_id"))
    }
  }
  
  attr(vals, "samples_dropped") <- dropped.samples
  attr(vals, "aggregated") <- aggregated
  set_fds(vals, x)
}


#' @noRd
#' @importFrom data.table as.data.table melt.data.table set setcolorder setnames
.melt.assay.matrix <- function(vals, assay_name, atype, ftype, finfo) {
  vals <- as.data.table(vals, keep.rownames=TRUE)
  vals <- melt.data.table(vals, id.vars='rn', variable.factor=FALSE,
                          variable.name='sample_id')
  setnames(vals, 1L, 'feature_id')
  
  set(vals, NULL, "dataset", sub("__.*$", "", vals[["sample_id"]]))
  set(vals, NULL, "sample_id", sub("^.*?__", "", vals[["sample_id"]]))
  set(vals, NULL, "assay", assay_name)
  set(vals, NULL, "assay_type", atype)
  set(vals, NULL, "feature_type", ftype)
  xref <- match(vals[["feature_id"]], finfo[["feature_id"]])
  set(vals, NULL, "feature_name", finfo[["name"]][xref])
  
  corder <- c("dataset", "sample_id", "assay", "assay_type",
              "feature_type", "feature_id", "feature_name", "value")
  setcolorder(vals, corder)
  vals
}

#' @noRd
#' @export
assay_names.FacileDataSet <- function(x, default_first = TRUE, ...) {
  anames <- assay_info_tbl(x) |> collect(n = Inf) |> pull(assay)
  if (default_first && length(anames) > 1L) {
    dassay <- default_assay(x)
    anames <- intersect(c(dassay, setdiff(anames, dassay)), anames)
  }
  anames
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
assay_sample_info.FacileDataSet <- function(x, assay_name, drop_samples = TRUE,
                                            samples = NULL, ...) {
  assert_facile_data_store(x)
  assert_choice(assay_name, assay_names(x))
  if (is.null(samples)) {
    samples <- collect(samples(x), n = Inf)
  } else {
    assert_sample_subset(samples, x)
  }
  
  asi <- assay_sample_info_tbl(x) |>
    filter(.data$assay == .env$assay_name) |>
    collect(n = Inf)
  
  dropped_samples <- NULL
  if (drop_samples) {
    dropped_samples <- anti_join(samples, asi, by = c("dataset", "sample_id"))
    out <- inner_join(samples, asi, by = c("dataset", "sample_id"), 
                      suffix = c(".x", ""))
  } else {
    out <- left_join(samples, asi, by = c("dataset", "sample_id"), 
                     suffix = c(".x", ""))  
  }
  
  attr(out, "samples_dropped") <- dropped_samples
  out
}
