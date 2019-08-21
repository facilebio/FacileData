#' Fetch data from single assay of choice
#'
#' @export
#' @importFrom rhdf5 h5read
#' @importFrom multiGSEA eigenWeightedMean
#' @rdname fetch_assay_data
#'
#' @md
#'
#' @param x A `FacileDataSrote` object.
#' @param features a feature descriptor (data.frame with assay and feature_id
#'   columms)
#' @param samples a sample descriptor to specify which samples to return data
#'   from.
#' @param assay_name the name of the assay to fetch data from. Defaults to the
#'   value of [default_assay()] for `x`. Must be a subset of `assay_names(x)`.
#' @param normalized return normalize or raw data values, defaults to
#'   \code{raw}
#' @param as.matrix by default, the data is returned in a long-form tbl-like
#'   result. If set to `TRUE`, the data is returned as a matrix.
#' @param ... parameters to pass to normalization methods
#' @param subset.threshold sometimes fetching all the genes is faster than
#'   trying to subset. We have to figure out why that is, but I've previously
#'   tested random features of different lengths, and around 700 features was
#'   the elbow.
#' @param aggregate.by do you want individual level results or geneset
#'   scores? Use 'ewm' for eigenWeightedMean, and that's all.
#' @return A `tibble` (lazy or not) with assay data.
#' @family API
fetch_assay_data.FacileDataSet <- function(x, features, samples = NULL,
                                           assay_name = default_assay(x),
                                           normalized = FALSE,
                                           as.matrix = FALSE,
                                           ...,
                                           subset.threshold = 700,
                                           aggregate = FALSE,
                                           aggregate.by= "ewm",
                                           verbose = FALSE) {
  assert_flag(as.matrix)
  assert_flag(normalized)
  assert_number(subset.threshold)

  if (!is.null(assay_name) || is.character(features)) {
    assert_string(assay_name)
    assert_choice(assay_name, assay_names(x))
  }

  if (missing(features) || is.null(features)) {
    assert_string(assay_name)
    features <- FacileData::features(x, assay_name) %>% collect(n=Inf)
  } else {
    if (is.character(features)) {
      features <- tibble(feature_id=features, assay=assay_name)
    }
    stopifnot(is(features, 'tbl') || is(features, 'data.frame'))
    if (!'assay' %in% colnames(features) || !is.character(features$assay)) {
      features <- collect(features, n = Inf)
      features[["assay"]] <- assay_name
    }
    assert_assay_feature_descriptor(features)
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

    assert_sample_subset(samples)
  }

  samples <- collect(samples)
  # you might think that you want to refactor this and put it before the
  # `collect()` call, but the SQLite back end can't handle
  # `distinct(..., .keep_all = TRUE)`
  samples <- distinct(samples, dataset, sample_id, .keep_all = TRUE)

  if (nrow(samples) == 0) {
    warning("Emtpy sample descriptor provided", immediate. = TRUE)
    return(samples)
  }

  assays <- unique(features$assay)
  n.assays <- length(assays)
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

  out <- lapply(assays, function(a) {
    f <- filter(features, assay == a)
    .fetch_assay_data(x, a, f$feature_id, samples, normalized, as.matrix,
                      subset.threshold, aggregate, aggregate.by, ...,
                      verbose = verbose)
  })

  if (length(out) == 1L) {
    out <- out[[1L]]
  } else {
    # We stop if we are asking for a matrix across multiple assays, but maybe
    # we don't have to ... (in the future, I mean)
    out <- bind_rows(out)
  }

  if (!as.matrix) {
    out <- as_facile_frame(out, x, classes = extra_classes,
                           .valid_sample_check = FALSE)
  }

  out
}

#' @export
#' @rdname fetch_assay_data
#' @importFrom data.table setDF
fetch_assay_data.facile_frame <- function(x, features, samples = NULL,
                                          assay_name = NULL,
                                          normalized = FALSE,
                                          as.matrix = FALSE,
                                          ...,
                                          subset.threshold = 700,
                                          aggregate = FALSE,
                                          aggregate.by= "ewm",
                                          verbose = FALSE) {
  fds. <- assert_facile_data_store(fds(x))
  if (!is.null(samples)) {
    warning("`samples` ignored when fetching covariates from a facile_frame",
            immediate. = TRUE)
  }
  samples. <- assert_sample_subset(x)
  if (is.null(assay_name)) assay_name <- default_assay(fds.)

  fetch_assay_data(fds., features = features, samples = samples.,
                   assay_name = assay_name, normalized = normalized,
                   as.matrix = as.matrix, ...,
                   subset.threshold = subset.threshold,
                   aggregate = aggregate, aggregate.by = aggregate.by,
                   verbose = verbose)
}


.fetch_assay_data <- function(x, assay_name, feature_ids, samples,
                              normalized = FALSE, as.matrix = FALSE,
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

  finfo <- features(x, assay_name, feature_ids = feature_ids) %>%
    collect(n = Inf) %>%
    arrange(hdf5_index)
  atype <- finfo$assay_type[1L]
  ftype <- finfo$feature_type[1L]
  sinfo <- assay_sample_info(x, assay_name, samples) %>%
    mutate(samid = paste(dataset, sample_id, sep = "__"))

  bad.samples <- is.na(sinfo$hdf5_index)
  if (any(bad.samples)) {
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

  dat <- sinfo %>%
    group_by(dataset) %>%
    do(res = {
      ds <- .$dataset[1L]
      hd5.name <- paste('assay', assay_name, ds, sep='/')
      vals <- h5read(hdf5fn(x), hd5.name, list(ridx, .$hdf5_index))
      if (is.null(ridx)) {
        vals <- vals[finfo$hdf5_index,,drop=FALSE]
      }
      dimnames(vals) <- list(finfo$feature_id, .$samid)
      if (normalized) {
        vals <- normalize.assay.matrix(vals, finfo, ., x, verbose=verbose, ...)
      }
      vals
    }) %>%
    ungroup

  # NOTE: We can avoid the monster matrix creation if we only want !as.matrix
  # returns, but this makes the code easier to reason about.
  # We can come back to this to optimize for speed later. The problem is
  # introduced when the aggregate.by parameter was introduced
  vals <- do.call(cbind, dat$res)

  if (nrow(vals) == 1L) {
    if (isTRUE(aggregate) && verbose) {
      warning("No assay feature aggregation performed over single feature",
              immediate.=TRUE)
    }
    aggregate <- FALSE
  }

  if (isTRUE(aggregate)) {
    scores <- switch(aggregate.by,
                     ewm = eigenWeightedMean(vals, ...)$score,
                     zscore = zScore(vals, ...)$score)
    vals <- matrix(scores, nrow=1, dimnames=list('score', names(scores)))
  }

  if (!as.matrix) {
    vals <- .melt.assay.matrix(vals, assay_name, atype, ftype, finfo)
    if (isTRUE(aggregate)) {
      vals[, feature_type := 'aggregated']
      vals[, feature_id := 'aggregated']
      vals[, feature_name := 'aggregated']
    }
    vals <- as.tbl(setDF(vals))
  }

  set_fds(vals, x)
}

#' @noRd
#' @importFrom data.table as.data.table melt.data.table set setcolorder setnames
.melt.assay.matrix <- function(vals, assay_name, atype, ftype, finfo) {
  vals <- as.data.table(vals, keep.rownames=TRUE)
  vals <- melt.data.table(vals, id.vars='rn', variable.factor=FALSE,
                          variable.name='sample_id')
  setnames(vals, 1L, 'feature_id')

  # vals[, dataset := sub('__.*$', '', sample_id)]
  # vals[, sample_id := sub('^.*__', '', sample_id)]
  # vals[, assay := assay_name]
  # vals[, assay_type := atype]
  # vals[, feature_type := ftype]
  # xref <- match(vals$feature_id, finfo$feature_id)
  # vals[, feature_name := finfo$name[xref]]

  set(vals, NULL, "dataset", sub("__.*$", "", vals[["sample_id"]]))
  set(vals, NULL, "sample_id", sub("^.*__", "", vals[["sample_id"]]))
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
    dat <- matrix(eigenWeightedMean(dat)$score, nrow = 1)
  }

}

#' @noRd
#' @export
assay_names.FacileDataSet <- function(x, default_first = TRUE, ...) {
  anames <- assay_info_tbl(x) %>% collect(n = Inf) %>% pull(assay)
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
  ainfo <- assay_info_tbl(x) %>% collect(n = Inf)
  if (!is.null(assay_name)) {
    assert_string(assay_name)
    assert_choice(assay_name, ainfo$assay)
    ainfo <- filter(ainfo, assay == assay_name)
  }
  as_facile_frame(ainfo, x)
}

#' @export
has_assay <- function(x, assay_name) {
  assert_facile_data_store(x)
  assert_character(assay_name)
  assay_name %in% assay_names(x)
}

#' Utility functions to get row and column indices of rnaseq hdf5 files.
#'
#' This is called to get things like hdf5_index and scaling factors for
#' the samples in a given assay.
#'
#' @export
#' @param x \code{FacileDataSet}
#' @param assay_name the name of the assay
#' @param samples a sample descriptor
#' @return an updated version of \code{samples} decorated with hd5_index,
#'   scaling factors, etc. Note that rows in \code{samples} that do not appear
#'   in \code{assay_name} will be returnd here with NA values for hd5_index and
#'   such.
assay_sample_info <- function(x, assay_name, samples = NULL) {
  assert_facile_data_store(x)

  if (is.null(samples)) {
    samples <- samples(x)
  } else {
    assert_sample_subset(samples, x)
    samples <- distinct(samples, dataset, sample_id, .keep_all = TRUE)
  }
  samples <- collect(samples, n = Inf)

  feature.type <- assay_feature_type(x, assay_name)

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


#' Returns the feature_type for a given assay
#'
#' The elements of the rows for a given assay all correspond to a particular
#' feature space (ie. feature_type='entrez')
#'
#' @export
#' @param x \code{FacileDataSet}
#' @param assay_name the name of the assay
assay_feature_type <- function(x, assay_name) {
  stopifnot(is.FacileDataSet(x))
  assert_string(assay_name)
  assert_choice(assay_name, assay_names(x))
  assay_info_tbl(x) %>%
    filter(assay == assay_name) %>%
    collect(n = Inf) %>%
    pull(feature_type)
}

#' Materializes a table with all feature information for a given assay.
#'
#' DEBUG: This logic is unnecessarily complex because I make sure to collect
#' all tables from the database as opposed to copying external tables in and
#' doing an inner_join in the database. I'm doing this becuase we are getting
#' name collisions on some of the temporary tables. We get errors like:
#'     Warning: Error in : Table pkdtpohpsu already exists.
#'
#' This fetches the hdf5_index for the assays as well
#'
#' @export
#' @inheritParams assay_feature_type
#' @param feature_ids a character vector of feature_ids
#' @return a \code{tbl_sqlite} result with the feature information for the
#'   features in a specified assay
assay_feature_info.FacileDataSet <- function(x, assay_name, feature_ids = NULL,
                                             ...) {
  .Deprecated("features()")
  features(x, assay_name, feature_ids, ...)
}

#' Retrieves feature information from the FacileDataStore for *either* a
#' particular `assay_name` or `feature_type`. By default this method returns
#' feature information for the features measured on the `default_assay()` of
#' this FacileDataStore.
#'
#' @export
#' @noRd
features.FacileDataSet <- function(x, assay_name = NULL, feature_type = NULL,
                                   feature_ids = NULL, ...) {
  null.aname <- is.null(assay_name)
  null.ftype <- is.null(feature_type)

  # Is the user asking for feature information from the features measured on
  # a given assay, or for all features of a given feature_type.
  if (null.aname && null.ftype) {
    assay_name <- default_assay(x)
    null.aname <- FALSE
  }
  if (!xor(null.aname, null.ftype)) {
    stop("Must specify feature information for EITHER assay_name or ",
         "feature_type")
  }
  if (null.aname) {
    query_type <- "feature_type"
    query_value <- assert_choice(feature_type, feature_types(x))
  } else {
    query_type <- "assay_name"
  }

  if (!is.null(feature_ids)) {
    assert_character(feature_ids)
  }

  if (query_type == "feature_type") {
    out <- filter(feature_info_tbl(x), feature_type == query_value)
    if (!is.null(feature_ids) && length(feature_ids) > 0) {
      out <- filter(out, feature_id %in% feautre_ids)
    }
    out <- collect(out, n = Inf)
  } else {
    ftype <- assay_feature_type(x, assay_name)

    afinfo <- assay_feature_info_tbl(x) %>%
      filter(assay == assay_name)

    if (!is.null(feature_ids) && length(feature_ids) > 0) {
      afinfo <- filter(afinfo, feature_id %in% feature_ids)
    }
    afinfo <- collect(afinfo, n=Inf)

    assay.info <- assay_info_tbl(x) %>%
      select(assay, assay_type, feature_type) %>%
      filter(assay == assay_name) %>%
      collect(n = Inf)

    ## FIXME: consider materialized view for this
    out <- inner_join(afinfo, assay.info, by = "assay")

    ftype <- out$feature_type[1L]
    finfo <- filter(feature_info_tbl(x), feature_type == ftype)
    finfo <- collect(finfo, n = Inf)

    # FIXME: feature_id should be made unique to feature_type to simplify
    # e.g add GeneID: prefix for entrez
    # But, still we know out and finfo each only have one feature type now
    out <- inner_join(out, finfo, by = c("feature_type", "feature_id"))
  }

  as_facile_frame(out, x)
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
#' @export
#' @param x FacileDataSet
#' @param samples sample descriptor
#' @param with_count return the number of samples in \code{samples} that are
#'   assayed over each assay as a column in \code{return}
#' @return rows from assay_info_tbl that correspond to the assays defined
#'   over the given samples. If no assays are defined over these samples,
#'   you're going to get an empty tibble.
assay_info_over_samples <- function(x, samples = NULL) {
  assert_facile_data_store(x)
  if (is.null(samples)) {
    samples <- samples(x)
  } else {
    assert_sample_subset(samples, x)
    samples <- distinct(samples, dataset, sample_id)
  }

  asi <- select(assay_sample_info_tbl(x), dataset, assay, sample_id)
  if (!same_src(asi, samples)) {
      asi <- collect(asi, n = Inf)
      samples <- collect(samples, n = Inf)
  }
  assays <- inner_join(asi, samples, by = c("dataset","sample_id"))

  # Count number of samples across dataset count for each assay type
  out <- assays %>%
    group_by(assay) %>%
    summarize(ndatasets = n_distinct(dataset), nsamples = n()) %>%
    ungroup() %>%
    collect(n = Inf)

  # default assay first
  def.idx <- which(out[["assay"]] == default_assay(x))
  if (length(def.idx)) {
    out <- bind_rows(out[def.idx,,drop = FALSE], out[-def.idx,,drop = FALSE])
  }

  as_facile_frame(out, x)
}

#' @section Removing Batch Effects:
#' We leverage limma's `removeBatchEffect` functionality, with a simplified
#' interface. The `batch` parameter replaces `batch`, `batch2`, and
#' `covariates`. The `design` parameter is replaced with `main`.
#'
#' All these parameters must be characters, which reference columns that
#' have been passed down into the samples frame, or ones that we can pull
#' out of the facile data store and spank onto the samples frame.
#'
#' We'll use these parameters to build a model.matrix with main and batch
#' effect and follow the use of `removeBatchEffect` as outlined in the post
#' linked to below to pull the design matrix apart and call the function with
#' the corresponding `design` and `covariates` parameters:
#'
#' https://support.bioconductor.org/p/83286/#83287
#'
#' Setting the `batch.scale` parameter to `TRUE` (the default), ensures that
#' the `rowMeans` of the returned data matrix are the same as the original
#' dataset.
#'
#' @rdname fetch_assay_data
#' @importFrom edgeR cpm
#' @importFrom limma removeBatchEffect
#' @importFrom stats contr.sum model.matrix
#' @examples
#' samples <- exampleFacileDataSet() %>%
#'   filter_samples(indication == "BLCA", sample_type == "tumor")
#' features <- c(PRF1='5551', GZMA='3001', CD274='29126')
#' dat <- with_assay_data(samples, features, normalized = TRUE, batch = "sex")
#' dat <- with_assay_data(samples, features, normalized = TRUE,
#'                        batch = c("sex", "RIN"), normalized = TRUE)
#' dat <- with_assay_data(samples, features, normalized = TRUE,
#'                        batch = c("sex", "RIN"), main = "treatment")
normalize.assay.matrix <- function(vals, feature.info, sample.info, fds,
                                   log = TRUE, prior.count = 1,
                                   batch = NULL, main = NULL,
                                   maintain.rowmeans = TRUE,
                                   verbose = FALSE, ...) {
  stopifnot(
    nrow(vals) == nrow(feature.info),
    all(rownames(vals) == feature.info$feature_id),
    ncol(vals) == nrow(sample.info),
    all(colnames(vals) == sample.info$samid),
    is.character(feature.info$assay_type),
    length(unique(feature.info$assay_type)) == 1L,
    is.numeric(sample.info$libsize), is.numeric(sample.info$normfactor))
  atype <- feature.info$assay_type[1L]
  if (atype %in% c("rnaseq", "isoseq")) {
    # we assume these are units that are at the count level
    # the user may have passes in a samples frame with lib.size and norm.factors
    # columns already attached, if so let's try using those
    # libsize <- sample.info$libsize * sample.info$normfactor
    cnames <- colnames(sample.info)
    if ("lib.size" %in% cnames && is.numeric(sample.info[["lib.size"]])) {
      lsize <- sample.info[["lib.size"]]
    } else  {
      lsize <- sample.info[["libsize"]]
    }
    if ("norm.factors" %in% cnames &&
        is.numeric(sample.info[["norm.factors"]])) {
      nf <- sample.info[["norm.factors"]]
    } else {
      nf <- sample.info[["normfactor"]]
    }
    libsize <- lsize * nf
    out <- edgeR::cpm(vals, libsize, log=log, prior.count=prior.count)
  } else if (atype == "tpm") {
    # someone processed their data with salmon or kallisto and wanted to store
    # tpm. Normalizing this is just log2(val + prior.count)
    out <- log2(vals + prior.count)
  } else {
    if (verbose) {
      warning("No normalization procedure for ", atype, " assay",
              immediate.=TRUE)
    }
    out <- vals
  }

  # We assume we've got log-transformed data here, and removeBatchEffect
  # if that was asked for.
  if (test_character(batch, min.len = 1L)) {
    if (is.character(main) && length(main) == 0L) main <- NULL
    if (!is.null(main)) {
      assert_string(main)
    }
    retrieve.covs <- setdiff(c(batch, main), colnames(sample.info))
    sample.info <- try({
      with_sample_covariates(sample.info, retrieve.covs, .fds = fds)
    }, silent = TRUE)
    if (is(sample.info, "try-error")) {
      stop("Covariates for batch correction could not be found: ",
           paste(retrieve.covs, collapse = ","))
    }

    # It's possible that the main and batch covariates are all singular,
    # ie. all factors with the same level, or whatever. Let's protect against
    # that before we removeBatchEffect
    batch.df <- sample.info[, c(main, batch), drop = FALSE]
    is.singular <- sapply(batch.df, function(vals) length(unique(vals)) == 1L)

    if (!is.null(main) && is.singular[main]) main <- NULL
    batch <- setdiff(batch, names(is.singular)[is.singular])

    if (length(batch)) {
      is.num <- sapply(sample.info[, batch, drop = FALSE], is.numeric)
      if (is.null(main)) {
        if (any(!is.num)) {
          cat.mats <- lapply(batch[!is.num], function(bcov) {
            # in limma we trust (code taken from limma::removeBatchEffect)
            batch. <- droplevels(as.factor(sample.info[[bcov]]))
            contrasts(batch.) <- contr.sum(levels(batch.))
            model.matrix(~ batch.)[, -1, drop = FALSE]
          })
          cat.mats <- do.call(cbind, cat.mats)
        } else {
          cat.mats <- matrix(0, nrow = nrow(sample.info), ncol = 0)
        }
        num.mats <- sample.info[, batch[is.num], drop = FALSE]
        batch.design <- cbind(cat.mats, num.mats)
        treatment.design <- matrix(1, nrow(sample.info), 1)
      } else {
        batch.formula <- paste(batch, collapse = " + ")
        des.formula <- paste("~", main, "+", batch.formula)
        des.matrix <- model.matrix(formula(des.formula), data = sample.info)
        main.cols <- c(1, grep(sprintf("^%s", main), colnames(des.matrix)))
        treatment.design <- des.matrix[, main.cols, drop = FALSE]
        batch.design <- des.matrix[, -(main.cols), drop = FALSE]
      }
      if (maintain.rowmeans) {
        batch.design <- scale(batch.design)
      }
      out <- removeBatchEffect(out, design = treatment.design,
                               covariates = batch.design)
    }
  }
  out
}

.extract_batch_terms <- function(sample.info, batch, batch2, covariates,
                                 design, batch_columns, ...) {
  # terms()
}

#' Creates a feature descriptor for interactive ease
#'
#' Creates a data.frame of features and assays they come from
#' @export
#' @param x FacileDataSet
#' @param features a character string of fearture ids (requires assay_name)
#'   or a data.frame with feature_id column.
#' @param assay_name the assay to get the featurespace from. if this is provided,
#'   it will trump an already existing assay_name column in \code{features}
#' @return a feature descriptor with feature_id and assay_name, which can be
#'   used to absolutely find features
create_assay_feature_descriptor <- function(x, features=NULL, assay_name=NULL) {
  ## TODO: Refactor the code inside `fetch_assay_data` to use this.
  # stopifnot(is.FacileDataSet(x))
  assert_facile_data_store(x)

  if (is.character(features) || is.null(features) || is(features, 'tbl_sql')) {
    if (is.null(assay_name)) assay_name <- default_assay(x)
    assert_string(assay_name)
    assert_choice(assay_name, assay_names(x))
  }

  if (is.null(features)) {
    features <- features(x, assay_name) %>% collect(n=Inf)
  } else if (is.character(features)) {
    features <- tibble(feature_id = features, assay = assay_name)
  } else if (is(features, 'tbl_sql')) {
    features <- mutate(collect(features, n = Inf), assay = assay_name)
  } else if (is.data.frame(features) && is.null(features[["assay"]])) {
    features[["assay"]] <- assay_name
  }

  assert_assay_feature_descriptor(features, x)
  features
}

#' Append expression values to sample-descriptor
#'
#' Since this is called in a "convenience" sort of way, often in a pipe-chain
#' \code{normalize} defaults to \code{TRUE}
#'
#' @export
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
  assert_flag(normalized)

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
  adata <- fetch_assay_data(.fds, features, x, normalized = normalized,
                            aggregate = aggregate, aggregate.by = aggregate.by,
                            verbose = verbose, ...)

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
    adata <- tidyr::spread_(adata, spread, 'value')
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

  out <- spread_(x, key, value) %>% set_fds(.fds)
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

  ss <- assay_sample_info_tbl(.fds) %>%
    filter(assay == assay_name) %>%
    select(dataset, sample_id, !!covariates) %>%
    collect(n = Inf)

  out <- left_join(x, ss, by = c("dataset", "sample_id"))
  as_facile_frame(out, ofds)
}
