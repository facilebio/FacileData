##' Converts a facile result into a DGEList or ExpressionSet, or ...
##'
##' The genes and samples that populate the \code{DGEList} are specified by
##' \code{x}, and the caller can request addition sample information to be
##' appended to \code{out$samples} via specification through the
##' \code{covariates} argument.
##'
##' @rdname expression-container
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
##' @rdname expression-container
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
##' @rdname expression-container
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
##' @rdname expression-container
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

##' @export
##' @rdname expression-container
as.ExpressionSet <- function(x, covariates=TRUE, feature_ids=NULL,
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

##' @export
##' @rdname expression-container
as.SummarizedExperiment <- function(x, covariates=TRUE, feature_ids=NULL,
                                    assay_name=default_assay(.fds),
                                    .fds=fds(x), custom_key=Sys.getenv("USER"),
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

## as.FacileDataSet conversion and utility functions ===========================

##' Converts bioconductor assay containers into a FacileDataSet.
##'
##' This function assumes you are only extracting one assay from the assay
##' container and creating a FacileDataSet from it. This requires that you
##' specify which assay (if the container has more than one) to extract, as
##' well as hand-crafting the feature_info correctly for the assay.
##'
##' @rdname faciledataset-converter
##' @export
##' @return a \code{\link{FacileDataSet}}
##' @param x The bioconductor assay container to extract data from
##' @param assay_name The name to use when storing the assay matrix from
##' \code{x} into the faciledataset.
##' @param assay_type what type of assay is this? rnaseq, microarry, nanostring,
##'   isoseq (isoform expression), etc.
##' @param feature_info a data.frame that describes the information for the
##'   features (rows) of the assay you are extracting. Currently you had to
##'   hand-craft this. In the future we will provide automated default
##'   fData, rowData, etc. extractors for the source assay cointaner.
##' @param feature_type \code{c('entrez', 'ensgid', 'enstid', 'genomic', 'custom')}
##' @param metayaml a yaml file (or list of lists) that describes the covariates
##'   in the pData/colData of \code{x}. If not provided, a default one will be
##'   generated
##' @param organism c("Homo sapiens", "Mus musculus", "unspecified"). This
##'   is used to fetch the appropriate genesets when this dataset is used with
##'   the facileexplorer
##' @param path the directory to create the faciledataset into. Will create
##'   a default directory in the current working directory if not specified.
##' @param source_assay the name of the assay element in \code{x} to extract
##'   for use.
##' @param ... more args
as.FacileDataSet <- function(x, assay_name, assay_type, feature_info,
                             feature_type, metayaml=NULL,
                             organism="unspecified",
                             path=tempfile("FacileDataSet-", getwd()),
                             source_assay=NULL, ...) {
  UseMethod("as.FacileDataSet")
}

##' @method as.FacileDataSet default
##' @export
as.FacileDataSet.default <- function(x, assay_name, assay_type, feature_info,
                                     feature_type, metayaml=NULL,
                                     organism="unspecified",
                                     path=tempfile("FacileDataSet-", getwd()),
                                     source_assay=NULL, ...) {
  stop("as.FacileDataSet not defined for object of class: ", class(x)[1L])
}

##' @method as.FacileDataSet ExpressionSet
##' @export
##' @rdname faciledataset-converter
as.FacileDataSet.ExpressionSet <- function(x, assay_name, assay_type, feature_info,
                                           feature_type, metayaml=NULL,
                                           organism="unspecified",
                                           path=tempfile("FacileDataSet-", getwd()),
                                           source_assay=assayDataElementNames(x)[1L],
                                           ...) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package required to convert ExpresionSet to FacileDataSet")
  }
}

##' @method as.FacileDataSet SummarizedExperiment
##' @export
##' @rdname faciledataset-converter
as.FacileDataSet.SummarizedExperiment <- function(x, assay_name, assay_type, feature_info,
                                                  feature_type, metayaml=NULL,
                                                  organism="unspecified",
                                                  path=tempfile("FacileDataSet-", getwd()),
                                                  source_assay=NULL, ...) {
}

##' @method as.FacileDataSet DGEList
##' @export
##' @rdname faciledataset-converter
as.FacileDataSet.DGEList <- function(x, assay_name, assay_type, feature_info,
                                     feature_type, metayaml=NULL,
                                     organism="unspecified",
                                     path=tempfile("FacileDataSet-", getwd()),
                                     source_assay=NULL, ...) {
}

as.FacileDataSet.matrix <- function(x, assay_name, assay_type, feature_info,
                                    feature_type, metayaml=NULL,
                                    organism="unspecified",
                                    path=tempfile("FacileDataSet-", getwd()),
                                    source_assay=NULL, ...) {

}

##' @method as.FacileDataSet list
##' @export
##' @rdname faciledataset-converter
as.FacileDataSet.list <- function(x, path, organism, assays=NULL, metayaml=NULL,
                                  ...) {
}

## Utlity functions to create meta.yaml file from various sample-covariates ====
metayaml_from_df <- function(x) {

}



#' Creates a shell of a yaml file for a FacileDataSet
create_metayaml <- function(name='unspecified', organism='unspecified',
                            datasets=list(), sample_covariates=list(),
                            default_assay='unspecified') {
  c('name', 'organism', 'datasets', 'sample_covariates',
    'default_assay')

}


metayaml_from_column <- function(x, columm, ...) {

}
