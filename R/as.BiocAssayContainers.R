##' Converts a facile result into a traditional Bioconductor assay containers.
##'
##' Even though the faciledata API provides functionality to access your data
##' in useful ways, it is still conceivable that you *might* prefer to work
##' with these data using a more traditional bioconductor container like a
##' `SummarizedExperiment`, `DGEList`, `ExpressionSet`, etc.These functions are
##' here for you.
##'
##' @rdname as.BiocContainer
##' @export
##' @importFrom edgeR DGEList
##' @param x a facile expression-like result
##' @param covariates The covariates the user wants to add to the $samples of
##'   the DGEList. This can take the following forms:
##'   \enumerate{
##'     \item{TRUE}{
##'       All covariates are retrieved from the \code{FacileDataSet}
##'     }
##'     \item{FALSE}{
##'       TODO: Better handle FALSE
##'     }
##'     \item{character}{
##'       A vector of covariate names to fetch from the \code{FacileDataSet}
##'     }
##'     \item{data.frame}{
##'       A table that looks like a subset of the sample_covariate table
##'     }
##'     \item{NULL}{
##'       If \code{NULL}, no covariates are retrieved
##'     }
##'   }
##' @param feature_ids the features to get expression for (if not specified
##'   in \code{x} descriptor)
##' @param assay the \code{assayDataElement} to use for the expression data
##'   if \code{x} is an \code{ExpressionSet}.
##' @param .fds The \code{FacileDataSet} that \code{x} was retrieved from.
##' @param custom_key the custom key to use to fetch custom annotations from
##'   \code{.fds}
##' @return a \code{\link[edgeR]{DGEList}}
as.DGEList <- function(x, ...) {
  UseMethod('as.DGEList')
}

##' @method as.DGEList matrix
##' @rdname as.BiocContainer
as.DGEList.matrix <- function(x, covariates=TRUE, feature_ids=NULL,
                              assay_name=default_assay(.fds), .fds=fds(x),
                              custom_key=Sys.getenv("USER"), ...) {
  ## NOTE: by now assay_name is ignored
  stopifnot(is(x, 'FacileExpression'))
  requireNamespace("edgeR")
  .fds <- force(.fds)
  stopifnot(is.FacileDataSet(.fds))

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
    } else if (is.character(covariates)) {
      covariates <- fetch_sample_covariates(.fds, samples, covariates)
    }
    assert_sample_covariates(covariates)
  }

  fids <- rownames(x)
  genes <- gene_info_tbl(.fds) %>%
    collect(n=Inf) %>% ## #dboptimize# remove this if you want to exercise db
    semi_join(tibble(feature_id=fids), by='feature_id') %>%
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
  sample.stats <- fetch_sample_statistics(.fds, samples) %>%
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
    covs <- spread_covariates(covariates, .fds) %>%
      as.data.frame %>%
      set_rownames(., paste(.$dataset, .$sample_id, sep='__')) %>%
      select(-dataset, -sample_id)
    y$samples <- cbind(y$samples, covs[colnames(y),,drop=FALSE])
  }

  set_fds(y, .fds)
}

##' @export
##' @method as.DGEList data.frame
##' @rdname as.BiocContainer
as.DGEList.data.frame <- function(x, covariates=TRUE, feature_ids=NULL,
                                  assay_name=default_assay(.fds), .fds=fds(x),
                                  custom_key=Sys.getenv("USER"),
                                  ...) {
  .fds <- force(.fds)
  stopifnot(is.FacileDataSet(.fds))
  x <- assert_sample_subset(x)

  has.count <- 'value' %in% colnames(x) && is.integer(x[['value']])
  fetch.counts <- !has.count

  ## Do we want to fetch counts from the FacileDataSet?
  if (has.count) {
    if (is.character(feature_ids) && all(feature_ids %in% x[['feature_id']])) {
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
    counts <- fetch_assay_data(.fds, feature_ids, x, assay_name=assay_name,
                               normalized=FALSE, as.matrix=TRUE)
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

  as.DGEList(counts, covariates=covariates, feature_ids=feature_ids,
             .fds=.fds, custom_key=custom_key, ...)
}

##' @export
##' @method as.DGEList tbl_sql
##' @rdname as.BiocContainer
as.DGEList.tbl_sql <- function(x, covariates=TRUE, feature_ids=NULL,
                               assay_name=default_assay(.fds), .fds=fds(x),
                               custom_key=Sys.getenv("USER"),
                               ...) {
  x <- collect(x, n=Inf) %>% set_fds(.fds)
  as.DGEList(x, covariates, feature_ids, assay_name, .fds=.fds,
             custom_key=custom_key, ...)
}

as.DGEList.FacileDataSet <- function(x, covariates=TRUE, feature_ids=NULL,
                                     assay_name=default_assay(x),
                                     custom_key=Sys.getenv("USER"),
                                     ...) {
  as.DGEList(samples(x), covariates, feature_ids, assay_name, x, custom_key,
             ...)
}

##' @rdname as.BiocContainer
##' @export
##' @return a \code{\link[Biobase]{ExpressionSet}}
as.ExpressionSet <- function(x, ...) {
  UseMethod('as.ExpressionSet')
}

##' @rdname as.BiocContainer
##' @export
##' @method as.ExpressionSet data.frame
##' @rdname as.BiocContainer
as.ExpressionSet.data.frame <- function(x, covariates=TRUE, feature_ids=NULL,
                                        assay_name=default_assay(.fds),
                                        .fds=fds(x), custom_key=Sys.getenv("USER"), ...) {
  .fds <- force(.fds)
  stopifnot(is.FacileDataSet(.fds))
  assert_sample_subset(x)
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase required")
  }
  y <- as.DGEList(x, covariates, feature_ids, assay_name, .fds=.fds,
                  custom_key=custom_key, ...)
  es <- Biobase::ExpressionSet(y$counts)
  es <- Biobase::`pData<-`(es, y$samples)
  es <- Biobase::`fData<-`(es, y$genes)
  set_fds(es, .fds)
}

##' @rdname as.BiocContainer
##' @export
##' @method as.ExpressionSet FacileDataSet
##' @rdname as.BiocContainer
as.ExpressionSet.FacileDataSet <- function(x, covariates=TRUE, feature_ids=NULL,
                                           assay_name=default_assay(.fds),
                                           .fds=fds(x),
                                           custom_key=Sys.getenv("USER"), ...) {
  force(.fds)
  x <- samples(x) %>% collect(n=Inf) %>% set_fds(.fds)
  as.ExpressionSet(x, covariates, feature_ids, assay_name, x,
                   custom_key, ...)
}

##' @rdname as.BiocContainer
##' @export
##' @return a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
as.SummarizedExperiment <- function(x, ...) {
  UseMethod('as.SummarizedExperiment')
}

##' @rdname as.BiocContainer
##' @export
##' @method as.SummarizedExperiment data.frame
##' @rdname as.BiocContainer
as.SummarizedExperiment.data.frame <- function(x, covariates=TRUE, feature_ids=NULL,
                                               assay_name=default_assay(.fds),
                                               .fds=fds(x),
                                               custom_key=Sys.getenv("USER"),
                                               ...) {
  .fds <- force(.fds)
  stopifnot(is.FacileDataSet(.fds))
  assert_sample_subset(x)
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package required")
  }
  y <- as.DGEList(x, covariates, feature_ids, assay_name, .fds=.fds,
                  custom_key=custom_key, ...)
  ## TODO: Check y$genes to see if we should make a rowRanges out of the
  ## rowData or just keep it as a DataFrame
  out <- SummarizedExperiment::SummarizedExperiment(
    y$counts, colData=y$samples, rowData=y$genes, ...)
  set_fds(out, .fds)
}

##' @rdname as.BiocContainer
##' @export
##' @method as.ExpressionSet FacileDataSet
##' @rdname as.BiocContainer
as.SummarizedExperiment.FacileDataSet <- function(x, covariates=TRUE, feature_ids=NULL,
                                                  assay_name=default_assay(.fds),
                                                  .fds=fds(x), custom_key=Sys.getenv("USER"),
                                                  ...) {
  force(.fds)
  x <- samples(x) %>% collect(n=Inf) %>% set_fds(.fds)
  as.SummarizedExperimentt(x, covariates, feature_ids, assay_name, x,
                           custom_key, ...)
}
