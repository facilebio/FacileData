#' Converts a "facile object" to a traditional Bioconductor assay container
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
#' @importFrom edgeR DGEList
as.DGEList.matrix <- function(x, covariates = TRUE, feature_ids = NULL,
                              assay_name = NULL, .fds = NULL,
                              custom_key = Sys.getenv("USER"), ...) {
  .fds <- assert_facile_data_store(.fds)
  assert_choice(assay_name, assay_names(.fds))

  ## Construct sample table from colnames of the matrix, and make sure this is
  ## legit
  samples <- tibble(
    dataset=sub('__.*$', '', colnames(x)),
    sample_id=sub('^.*?__', '', colnames(x)))
  ## if you don't want to `collect` first, you could send `samples` in as
  ## second argument and then copy that into the db.
  ## #dboptimize
  bad.samples <- samples %>%
    anti_join(collect(sample_stats_tbl(.fds), n=Inf),
              by=c('dataset', 'sample_id')) %>%
    collect(n=Inf)
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

  # Construct $genes meta information
  ainfo <- assay_info(.fds, assay_name = assay_name)
  ftype <- ainfo[["feature_type"]]
  fids <- rownames(x)

  gene_info <- feature_info_tbl(.fds) %>%
    filter(feature_type == !!ftype) %>%
    collect(n = Inf)
  genes <- gene_info %>%
    semi_join(tibble(feature_id=fids), by='feature_id') %>%
    rename(symbol = "name") %>%
    as.data.frame %>%
    set_rownames(., .$feature_id)

  class(x) <- 'matrix'

  ## now subset down to only features asked for
  if (!is.null(feature_ids) && is.character(feature_ids)) {
    keep <- feature_ids %in% rownames(x)
    if (mean(keep) != 1) {
      warning(sprintf("Only %d / %d feature_ids requested are in dataset",
                      sum(keep), length(keep)))
    }
    x <- x[feature_ids[keep],,drop=FALSE]
    genes <- genes[feature_ids[keep],,drop=FALSE]
  }

  ## Doing the internal filtering seems to be too slow
  ## sample.stats <- fetch_sample_statistics(db, x) %>%
  sample.stats <- .fds %>%
    fetch_sample_statistics(samples, assay_name = assay_name) %>%
    collect(n=Inf) %>%
    mutate(samid=paste(dataset, sample_id, sep='__')) %>%
    rename(lib.size=libsize, norm.factors=normfactor) %>%
    as.data.frame %>%
    set_rownames(., .$samid)
  sample.stats <- sample.stats[colnames(x),,drop=FALSE]

  y <- DGEList(x, genes=genes, lib.size=sample.stats$lib.size,
               norm.factors=sample.stats$norm.factors)

  y$samples <- cbind(
    y$samples,
    sample.stats[colnames(y), c('dataset', 'sample_id', 'samid'), drop=FALSE])

  if (!is.null(covariates)) {
    covs <- select(covariates, -dataset, -sample_id)
    y$samples <- cbind(y$samples, covs[colnames(y),,drop=FALSE])
  }

  set_fds(y, .fds)
}

#' @export
#' @method as.DGEList data.frame
#' @rdname as.BiocContainer
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
    if (ainfo$assay_type != 'rnaseq') {
      warning("Creating DGEList for something other than rnaseq type assay")
    }
    counts <- fetch_assay_data(.fds, feature_ids, x, assay_name = assay_name,
                               normalized = FALSE, as.matrix = TRUE)
  } else {
    counts.dt <- assert_expression_result(x) %>%
      collect(n=Inf) %>%
      setDT %>%
      unique(by=c('dataset', 'sample_id', 'feature_id'))
    counts.dt[, samid := paste(dataset, sample_id, sep='__')]
    counts <- local({
      wide <- dcast.data.table(counts.dt, feature_id ~ samid, value.var='value')
      out <- as.matrix(wide[, -1L, with=FALSE])
      rownames(out) <- wide[[1L]]
      class(out) <- c('FacileExpression', class(out))
      out
    })
  }

  as.DGEList.matrix(counts, covariates = covariates, feature_ids = feature_ids,
                    assay_name = assay_name, .fds = .fds,
                    custom_key = custom_key, ...)
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

  # force(.fds)
  # force(assay_name)
  # .fds <- assert_facile_data_store(.fds)
  # browser()
  # NextMethod(.fds = .fds)
  # NextMethod(.fds = .fds)
  # NextMethod()
  # browser()

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
