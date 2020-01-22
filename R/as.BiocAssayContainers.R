#' Converts a "facile result" to a traditional Bioconductor assay container.
#'
#' An entire `FacileDataSet` or a subset of it can be converted into
#' bioconductor-standard assay containers, like a `SummarizedExperiment`,
#' `DGEList`, or `ExpressionSet` "at any time" using various `as.XXX` functions,
#' like `as.DGEList(...)`.
#'
#' We use the term "facile object" to refer to either the entirety of a
#' `FacileDataStore` or any sample-descriptor that specifies subsets of the
#' data, eg. where `fds(x)` returns a `FacileDataStore`. See examples for
#' specifics.
#'
#' Note that the order that the samples and features are materialized into the
#' expression container are not guaranteed.
#'
#' @rdname as.BiocContainer
#'
#' @export
#' @importFrom edgeR DGEList
#'
#' @param x a facile expression-like result
#' @param covariates The covariates the user wants to add to the $samples of
#'   the DGEList. This can take the following forms:
#'   - `TRUE`: All covariates are retrieved from the `FacileDataSet`
#'   - `FALSE`: TODO: Better handle FALSE
#'   - `character`: A vector of covariate names to fetch from the
#'     `FacileDataSet`. Must be elements of `names(sample_definitions(x))`
#'   - `data.frame`: A wide covariate table (dataset, sample_id, covariates ...)
#'     This may be external covariates for samples not available within
#'     `x` (yet), ie. a table of covariates provided by a third party.
#'   - `NULL`: do not decorate with *any* covariates.
#' @param feature_ids the features to get expression for (if not specified
#'   in `x` descriptor). These correspond to the elements found in the
#'   `feature_info_tbl(x)$feature_id` column.
#' @param assay_name the name of the assay matrix to use when populating the
#'   default assay matrix of the bioconductor container (the `$counts`
#'   matrix of a `DGEList`, the `exprs()` of an `ExpressionSet`, etc.).
#'   The default value is the entry provided by [default_assay()]
#' @param .fds The `FacileDataSet` that `x` was retrieved from
#' @param custom_key the custom key to use to fetch custom annotations from
#'   `.fds`
#' @return the appropriate bioconductor assay container, ie. an `edgeR::DGEList`
#'   for `as.DGEList`, a `Biobase::ExpressionSet` for `as.ExpressionSet`, or
#'   a `SummarizedExperiment::SummarizedExperiment` for
#'   `as.SummarizedExperiment`.
#'
#' @examples
#' fds <- exampleFacileDataSet()
#'
#' # Retrieve DGEList of gene expression for all samples
#' y.all <- as.DGEList(fds) # gene expression of all samples
#'
#' # Retrieve data for only 3 genes
#' # Suppose we only wanted female samples in our DGEList
#' y.fem <- fds %>%
#'   filter_samples(sex == "f") %>%
#'   as.DGEList() # or `as.ExpressionSet()`
#' @export
as.DGEList <- function(x, ...) {
  UseMethod("as.DGEList", x)
}

#' @noRd
#' @method as.DGEList matrix
#' @rdname as.BiocContainer
#' @importFrom edgeR DGEList calcNormFactors
as.DGEList.matrix <- function(x, covariates = TRUE, feature_ids = NULL,
                              assay_name = NULL, .fds = NULL,
                              custom_key = Sys.getenv("USER"),
                              update_libsizes = !is.null(feature_ids),
                              update_normfactors = update_libsizes,
                              custom_libsizes = FALSE,
                              custom_normfactors = FALSE,
                              ...) {
  .fds <- assert_facile_data_store(.fds)
  assert_choice(assay_name, assay_names(.fds))

  ## Construct sample table from colnames of the matrix, and make sure this is
  ## legit
  samples <- tibble(
    dataset=sub('__.*$', '', colnames(x)),
    sample_id=sub('^.*?__', '', colnames(x)))

  # if you don't want to `collect` first, you could send `samples` in as
  # second argument and then copy that into the db.
  # #dboptimize

  # TODO: Start fixing here. This msethod assumes direct access to the db of
  # the .fds, which breaks the abstractions. Fix use of these methods:
  #   * feature_info_tbl()
  #   * fetch_sample_statistics() should be fetch_assay_covariates? or
  #     with_assay_covariates, because we get a wide table, maybe adding
  #     the libsize and normfactors would look like:
  #     assay_info <- with_assay_covariates(samples, assay_name)
  bad.samples <- samples %>%
    anti_join(collect(samples(.fds), n = Inf),
              by=c("dataset", "sample_id")) %>%
    collect(n = Inf)
  if (nrow(bad.samples)) {
    stop("Bad sample columns specified in the count matrix")
  }

  ## Fetch appropriate covariate
  if (!is.null(covariates)) {
    if (isTRUE(covariates)) {
      covariates <- fetch_sample_covariates(.fds, samples)
      covariates <- spread_covariates(covariates)
    } else if (is.character(covariates)) {
      covariates <- fetch_sample_covariates(.fds, samples, covariates)
      covariates <- spread_covariates(covariates)
    }
    assert_subset(c("dataset", "sample_id"), colnames(covariates))
    assert_true(nrow(covariates) == nrow(distinct(covariates)))
    covariates <- as.data.frame(covariates, stringsAsFactors = FALSE)
    rownames(covariates) <-  paste(covariates$dataset,
                                   covariates$sample_id,
                                   sep="__")
  }

  if (is(feature_ids, "data.frame") || is(feature_ids, "tbl")) {
    feature_ids <- feature_ids[["feature_id"]]
  }
  if (!is.null(feature_ids) && is.character(feature_ids)) {
    keep <- feature_ids %in% rownames(x)
    if (mean(keep) != 1) {
      warning(sprintf("Only %d / %d feature_ids requested are in data matrix",
                      sum(keep), length(keep)))
    }
    feature_ids <- feature_ids[keep]
  } else {
    feature_ids <- rownames(x)
  }

  genes <- .fds %>%
    features(assay_name = assay_name, feature_ids = feature_ids) %>%
    as.data.frame()
  rownames(genes) <- genes[["feature_id"]]
  if (genes[["feature_type"]][1L] %in% c("entrez", "ensgid") || TRUE) {
    # After adding fluidigm/qPCR support, I wish I never renamed `name` to
    # `symbol`
    genes <- rename(genes, symbol = "name")
  }

  class(x) <- 'matrix'
  x <- x[feature_ids,,drop = FALSE]
  genes <- genes[feature_ids,,drop = FALSE]

  # Issue #2
  # https://github.com/denalitherapeutics/FacileData/issues/2
  sample.stats <- samples %>%
    with_assay_covariates(c("libsize", "normfactor"), assay_name,
                          .fds = .fds) %>%
    collect(n = Inf) %>%
    mutate(samid=paste(dataset, sample_id, sep='__')) %>%
    as.data.frame()
  rownames(sample.stats) <- sample.stats[["samid"]]
  sample.stats <- sample.stats[colnames(x),,drop=FALSE]

  # TODO: https://github.com/facileverse/FacileData/issues/9
  # Since this is being used a general assay retrieval wrapper, what we
  # want might not even be appropriate for a DGEList (ie. negative numbers
  # inside) we hack for now, but need to address this!
  is.neg <- which(x < 0, arr.ind = TRUE)
  if (nrow(is.neg)) x[is.neg] <- x[is.neg] * -1

  y <- DGEList(x, genes = genes, lib.size = sample.stats[["libsize"]],
               norm.factors = sample.stats[["normfactor"]])

  y$samples <- cbind(
    y$samples,
    sample.stats[, c("dataset", "sample_id", "samid"), drop = FALSE])

  if (update_libsizes) {
    y[["samples"]][["lib.size"]] <- colSums(y[["counts"]])
  }
  if (update_normfactors) {
    y <- edgeR::calcNormFactors(y)
  }

  if (!is.null(covariates)) {
    covs <- select(covariates, -dataset, -sample_id)
    covs <- covs[colnames(y),,drop=FALSE]
    if ("group" %in% colnames(covs)) {
      y[["samples"]][["group"]] <- covs[["group"]]
      covs[["group"]] <- NULL
    }
    # Remove any lib.size or norm.factors that might have come in through
    # here. Replace the ones loaded in the DGEList if the corresponding
    # `custom_*` parameter is set to `TRUE`
    if ("lib.size" %in% colnames(covs)) {
      if (custom_libsizes && is.numeric(covs[["lib.size"]])) {
        y$samples[["lib.size"]] <- covs[["lib.size"]]
      }
      covs[["lib.size"]] <- NULL
    }
    if ("norm.factors" %in% colnames(covs)) {
      if (custom_normfactors && is.numeric(covs[["norm.factors"]])) {
        y$samples[["norm.factors"]] <- covs[["norm.factors"]]
      }
      covs[["norm.factors"]] <- NULL
    }

    y$samples <- cbind(y$samples, covs)
  }

  # TODO: https://github.com/facileverse/FacileData/issues/9
  if (nrow(is.neg)) {
    y$counts[is.neg] <- y$counts[is.neg] * -1
  }

  set_fds(y, .fds)
}

#' @export
#' @method as.DGEList data.frame
#' @rdname as.BiocContainer
#' @importFrom data.table dcast set setDT
as.DGEList.data.frame <- function(x, covariates = TRUE, feature_ids = NULL,
                                  assay_name = NULL, .fds = NULL,
                                  custom_key = Sys.getenv("USER"),
                                  ...) {
  .fds <- assert_facile_data_store(.fds)
  if (is.null(assay_name)) {
    assay_name <- default_assay(.fds)
  }
  assert_choice(assay_name, assay_names(.fds))
  x <- assert_sample_subset(x)

  has.count <- 'value' %in% colnames(x) && is.integer(x[['value']])
  fetch.counts <- !has.count

  ## Do we want to fetch counts from the FacileDataSet?
  if (has.count) {
    if (is.character(feature_ids) && !all(feature_ids %in% x[['feature_id']])) {
      fetch.counts <- TRUE
    }
    if (!missing(feature_ids) && is.null(feature_ids)) {
      ## user explicitly wants everythin
      fetch.counts <- TRUE
    }
  }

  if (fetch.counts) {
    if (has.count) {
      warning("Ignoring expression in `x` and fetching data for `feature_ids`",
              immediate.=TRUE)
    }
    ## Check that we are getting the right type of assay for this
    ainfo <- assay_info(.fds, assay_name)
    if (!ainfo$assay_type %in% c("rnaseq", "isoseq")) {
      warning("Creating DGEList for something other than rnaseq type assay")
    }
    counts <- fetch_assay_data(.fds, feature_ids, x, assay_name = assay_name,
                               normalized = FALSE, as.matrix = TRUE)
  } else {
    counts.dt <- assert_expression_result(x) %>%
      collect(n = Inf) %>%
      setDT %>%
      unique(by = c('dataset', 'sample_id', 'feature_id'))
    # counts.dt[, samid := paste(dataset, sample_id, sep='__')]
    set(counts.dt, i = NULL, j = "samid",
        paste(counts.dt[["dataset"]], counts.dt[["sample_id"]], sep = "__"))
    counts <- local({
      wide <- dcast(counts.dt, feature_id ~ samid, value.var = "value")
      out <- as.matrix(wide[, -1L, with=FALSE])
      rownames(out) <- wide[[1L]]
      class(out) <- c('FacileExpression', class(out))
      out
    })
  }

  out <- as.DGEList.matrix(counts, covariates = covariates,
                           feature_ids = feature_ids,
                           assay_name = assay_name, .fds = .fds,
                           custom_key = custom_key, ...)

  # Add any covariates in x back to out$samples
  xref <- match(colnames(out), paste(x$dataset, x$sample_id, sep = "__"))
  cnames <- setdiff(colnames(x),
                    c("dataset", "sample_id", "lib.size", "norm.factors"))
  if (length(cnames)) {
    for (cname in cnames) out$samples[[cname]] <- x[[cname]][xref]
  }
  out
}

#' @method as.DGEList tbl
#' @export
#' @rdname as.BiocContainer
as.DGEList.tbl <- function(x, covariates = TRUE, feature_ids = NULL,
                           assay_name = NULL, .fds = NULL,
                           custom_key = Sys.getenv("USER"),
                           ...) {
  .fds <- assert_facile_data_store(.fds)
  if (is.null(assay_name)) {
    assay_name <- default_assay(.fds)
  }
  assert_choice(assay_name, assay_names(.fds))
  x <- collect(x, n = Inf)
  # NextMethod()
  as.DGEList.data.frame(x, covariates, feature_ids, assay_name, .fds = .fds,
                        custom_key = custom_key, ...)
}

#' @export
#' @rdname as.BiocContainer
as.DGEList.facile_frame <- function(x, covariates = TRUE, feature_ids = NULL,
                                    assay_name = NULL,
                                    custom_key = Sys.getenv("USER"),
                                    ...) {
  x <- collect(x, n = Inf)
  .fds <- assert_facile_data_store(fds(x))
  if (is.null(assay_name)) {
    assay_name <- default_assay(.fds)
  }
  assert_choice(assay_name, assay_names(.fds))

  has.count <- "value" %in% colnames(x) &&
    is.integer(x[["value"]]) &&
    is.character("feature_id")

  if (has.count && is.null(feature_ids)) {
    feature_ids <- unique(x[["feature_id"]])
  }

  as.DGEList.tbl(x, covariates, feature_ids, assay_name, .fds = .fds,
                 custom_key = custom_key, ...)
}


#' @export
#' @rdname as.BiocContainer
as.DGEList.FacileDataSet <- function(x, covariates = TRUE, feature_ids = NULL,
                                     assay_name = NULL,
                                     custom_key = Sys.getenv("USER"),
                                     ...) {
  xs <- samples(x)
  if (is.null(assay_name)) {
    assay_name <- default_assay(x)
  }
  assert_choice(assay_name, assay_names(x))
  as.DGEList(xs, covariates = covariates, feature_ids = feature_ids,
             assay_name = assay_name, custom_key = custom_key, ...)
}

#' @rdname as.BiocContainer
#' @export
#' @return a \code{\link[Biobase]{ExpressionSet}}
as.ExpressionSet <- function(x, ...) {
  UseMethod('as.ExpressionSet')
}

#' @rdname as.BiocContainer
#' @export
#' @method as.ExpressionSet data.frame
#' @rdname as.BiocContainer
as.ExpressionSet.data.frame <- function(x, covariates=TRUE, feature_ids=NULL,
                                        assay_name=default_assay(.fds),
                                        .fds=fds(x), custom_key=Sys.getenv("USER"), ...) {
  .fds <- force(.fds)
  stopifnot(is.FacileDataSet(.fds))
  assert_sample_subset(x)

  ns <- tryCatch(loadNamespace("Biobase"), error = function(e) NULL)
  if (is.null(ns)) stop("Biobase required for `as.ExpressionSet`")

  y <- as.DGEList(x, covariates, feature_ids, assay_name, .fds=.fds,
                  custom_key=custom_key, ...)
  es <- ns$ExpressionSet(y$counts)
  es <- ns$`pData<-`(es, y$samples)
  es <- ns$`fData<-`(es, y$genes)
  set_fds(es, .fds)
}

#' @rdname as.BiocContainer
#' @export
#' @method as.ExpressionSet FacileDataSet
#' @rdname as.BiocContainer
as.ExpressionSet.FacileDataSet <- function(x, covariates=TRUE, feature_ids=NULL,
                                           assay_name=default_assay(.fds),
                                           .fds=fds(x),
                                           custom_key=Sys.getenv("USER"), ...) {
  force(.fds)
  x <- samples(x) %>% collect(n=Inf) %>% set_fds(.fds)
  as.ExpressionSet(x, covariates, feature_ids, assay_name, x,
                   custom_key, ...)
}

#' @rdname as.BiocContainer
#' @export
#' @return a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
as.SummarizedExperiment <- function(x, ...) {
  UseMethod('as.SummarizedExperiment')
}

#' @rdname as.BiocContainer
#' @export
#' @method as.SummarizedExperiment data.frame
#' @rdname as.BiocContainer
as.SummarizedExperiment.data.frame <- function(x, covariates=TRUE, feature_ids=NULL,
                                               assay_name=default_assay(.fds),
                                               .fds=fds(x),
                                               custom_key=Sys.getenv("USER"),
                                               ...) {
  .fds <- force(.fds)
  stopifnot(is.FacileDataSet(.fds))
  assert_sample_subset(x)

  ns <- tryCatch(loadNamespace("SummarizedExperiment"), error = function(e) NULL)
  if (is.null(ns)) stop("SummarizedExperiment required for")

  y <- as.DGEList(x, covariates, feature_ids, assay_name, .fds=.fds,
                  custom_key=custom_key, ...)
  ## TODO: Check y$genes to see if we should make a rowRanges out of the
  ## rowData or just keep it as a DataFrame
  out <- ns$SummarizedExperiment(
    y$counts, colData=y$samples, rowData=y$genes, ...)
  set_fds(out, .fds)
}

#' @rdname as.BiocContainer
#' @export
#' @method as.SummarizedExperiment FacileDataSet
#' @rdname as.BiocContainer
as.SummarizedExperiment.FacileDataSet <- function(x, covariates=TRUE, feature_ids=NULL,
                                                  assay_name=default_assay(.fds),
                                                  .fds=fds(x), custom_key=Sys.getenv("USER"),
                                                  ...) {
  force(.fds)
  x <- samples(x) %>% collect(n=Inf) %>% set_fds(.fds)
  as.SummarizedExperiment(x, covariates, feature_ids, assay_name, x,
                           custom_key, ...)
}
