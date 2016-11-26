##' Utility functions to get row and column indices of rnaseq hdf5 files.
##'
##' @export
##' @rdname hdf5expression
hdf5_sample_indices <- function(x, samples=NULL) {
  stopifnot(is.FacileDataSet(x))
  if (is.null(samples)) {
    samples <- sample_stats_tbl(x)
  }
  samples <- samples %>%
    assert_sample_subset %>%
    collect(n=Inf) %>%
    distinct(dataset, sample_id)

  idxs <- hdf5_sample_xref_tbl(x) %>%
    collect(n=Inf) %>%
    semi_join(samples, by=c('dataset', 'sample_id')) %>%
    rename(index=hdf5_index)
}

##' @export
##' @rdname hdf5expression
hdf5_gene_indices <- function(x, feature_ids=NULL) {
  stopifnot(is.FacileDataSet(x))
  if (is.null(feature_ids)) {
    feature_ids <- character()
  }
  assertCharacter(feature_ids)
  feature_ids <- unique(feature_ids)

  genes <- gene_info_tbl(x)
  if (length(feature_ids) == 1L) {
    genes <- filter(genes, feature_id == feature_ids)
  } else if (length(feature_ids) > 1L) {
    genes <- filter(genes, feature_id %in% feature_ids)
  }

  genes %>%
    collect(n=Inf) %>%
    rename(index=hdf5_index)
}

##' Helper function creates dplyr query to get expression data of interest.
##'
##' @section Software Design Caveats:
##'
##' Note that if \code{db} is not provided, the result will be \code{collect}ed,
##' because the database connection will be closed when the function exits.
##'
##' If a \code{samples} selector is provided, we currently do not make sure that
##' each "observation" (dataset,sample_id)-tuple is unique. Not sure if
##' duplicated gene expression results will matter downstream ...
##'
##' @export
##' @importFrom rhdf5 h5read
##' @param x A \code{FacileDataSet} object.
##' @param indication The indication to get data from
##' @param subtype The subtype from the desired \code{indication}
##' @param sample_ids The sample_id's (barcodes) to fetch. If specified
##'   this trumps whatever is specified in \code{indication} and
##'   \code{subtype}
##' @param entrez_ids Restricts the fetch to only these genes
##' @return A lazy \code{\link[dplyr]{tbl}} object with the expression
##'   data to be \code{\link[dplyr]{collect}}ed when \code{db} is provided,
##'   othwerise a \code{tbl_df} of the results.
fetch_expression <- function(x, samples=NULL, feature_ids=NULL,
                             with_symbols=!as.matrix, as.matrix=FALSE) {
  stopifnot(is.FacileDataSet(x))

  hdf.sample.idxs <- hdf5_sample_indices(x, samples)
  hdf.gene.idxs <- hdf5_gene_indices(x, feature_ids)

  if (isTRUE(as.matrix) && isTRUE(with_symbols)) {
    warning("with_symbols ignored when as.matrix=TRUE")
  }

  counts <- hdf.sample.idxs %>%
    group_by(dataset) %>%
    do(res={
      ds <- .$dataset[1L]
      hd5.name <- paste0('expression/rnaseq/', ds)
      cnts <- h5read(x$hdf5.fn, hd5.name, list(hdf.gene.idxs$index, .$index))
      if (isTRUE(as.matrix)) {
        dimnames(cnts) <- list(hdf.gene.idxs$feature_id,
                               paste(ds, .$sample_id, sep='_'))
      } else {
        dimnames(cnts) <- list(hdf.gene.idxs$feature_id, .$sample_id)
        cnts <- as.data.table(cnts, keep.rownames=TRUE)
        cnts <- melt.data.table(cnts, id.vars='rn')
        cnts[, dataset := ds]
        cnts[, variable := as.character(variable)]
        setnames(cnts, c('feature_id', 'sample_id', 'count', 'dataset'))
        setcolorder(cnts, c('dataset', 'sample_id', 'feature_id', 'count'))
        if (isTRUE(with_symbols)) {
          sxref <- match(cnts$feature_id, hdf.gene.idxs$feature_id)
          cnts[, symbol := hdf.gene.idxs$symbol[sxref]]
        }
        cnts
      }
      cnts
    }) %>%
    ungroup

  if (isTRUE(as.matrix)) {
    out <- do.call(cbind, counts$res)
  } else {
    out <- as.tbl(setDF(rbindlist(counts$res)))
  }

  class(out) <- c('FacileExpression', class(out))
  set_fds(out, x)
}

##' Append expression values to sample-descriptor
##'
##' @export
##' @param x a samples descriptor
##' @param feature_ids character vector of feature_ids
##' @param with_symbols Do you want gene symbols returned, too?
##' @param .fds A \code{FacileDataSet} object
##' @return a tbl-like result
with_expression <- function(samples, feature_ids, with_symbols=TRUE,
                            .fds=fds(samples)) {
  stopifnot(is.FacileDataSet(.fds))
  stopifnot(is.character(feature_ids) && length(feature_ids) > 0)
  samples <- assert_sample_subset(samples)

  .fds %>%
    fetch_expression(samples, feature_ids, with_symbols=with_symbols) %>%
    join_samples(samples) %>%
    set_fds(.fds)
}

##' @importFrom edgeR cpm
##' @export
cpm <- function(x, ...) UseMethod("cpm")

##' Calculated counts per million from a FacileDataSet
##'
##' @method cpm tbl_sqlite
##' @rdname cpm
##'
##' @importFrom edgeR cpm
##' @export
##'
##' @param x an expression-like facile result
##' @param lib.size ignored for now, this is fetched from the
##'   \code{FacileDataSet}
##' @param log log the result?
##' @param prior.count prior.count to add to observed counts
##' @param .fds A \code{FacileDataSet} object
##' @return a modified expression-like result with a \code{cpm} column.
cpm.tbl_sqlite <- function(x, lib.size=NULL, log=FALSE, prior.count=5,
                           .fds=fds(x), ...) {
  assert_expression_result(x)
  stopifnot(is.FacileDataSet(.fds))
  if (!is.null(lib.size)) {
    warning("not supporting custom lib.size yet")
  }
  ## Let's leverage the fact that we're already working in the database
  sample.stats <- fetch_sample_statistics(.fds)

  if (same_src(x, sample.stats)) {
    tmp <- inner_join(x, sample.stats, by=c('dataset', 'sample_id'))
    tmp <- collect(tmp, n=Inf)
    cpms <- calc.cpm(tmp, lib.size=tmp$libsize * tmp$normfactor,
                     log=log, prior.count=prior.count)
    out <- tmp %>%
      mutate(cpm=cpms) %>%
      select_(.dots=c(colnames(x), 'cpm'))
  } else {

    out <- cpm(collect(x, n=Inf), lib.size=lib.size, log=log,
               prior.count=prior.count, sample.stats=sample.stats, .fds=.fds,
               ...)

  }
  set_fds(out, .fds)
}

##' @method cpm tbl_df
##' @rdname cpm
##' @export
cpm.tbl_df <- function(x, lib.size=NULL, log=FALSE, prior.count=5,
                       sample.stats=NULL, .fds=fds(x), ...) {
  assert_expression_result(x)
  stopifnot(is.FacileDataSet(.fds))
  if (!is.null(lib.size)) {
    warning("not supporting custom lib.size yet")
  }

  if (is.null(sample.stats)) {
    samples <- distinct(x, dataset, sample_id)
    sample.stats <- fetch_sample_statistics(.fds, samples)
  }

  sample.stats <- sample.stats %>%
    assert_sample_statistics %>%
    collect(n=Inf) %>%
    distinct(dataset, sample_id, .keep_all=TRUE)
  # if (is.data.table(x)) {
  #   setDT(sample.stats)
  # }

  tmp <- inner_join(x, sample.stats, by=c('dataset', 'sample_id'))
  cpms <- calc.cpm(tmp, lib.size=tmp$libsize * tmp$normfactor, log=log,
                   prior.count=prior.count)

  tmp %>%
    mutate(cpm=cpms) %>%
    select_(.dots=c(colnames(x), 'cpm')) %>%
    set_fds(.fds)
}


## @method cpm tbl_dt
## @export
# cpm.tbl_dt <- function(x, lib.size=NULL, log=FALSE, prior.count=5,
#                        sample.stats=NULL, .fds=fds(x), ...) {
#   cpm.tbl_df(x, lib.size, log, prior.count, sample.stats, .fds, ...)
# }


## A helper function to calculate cpm. Expects that x is a tbl-like thing which
## has a count column, and the appropriate libsize and normfactor columns from
## the database
calc.cpm <- function(x, lib.size=NULL, log=FALSE, prior.count=5) {
  x <- collect(x, n=Inf)

  pr.count <- lib.size / mean(lib.size) * prior.count

  if (log) {
    lib.size <- 1e-6 * (lib.size + 2*pr.count)
    cpms <- log2((x$count + pr.count) / lib.size)
  } else {
    lib.size <- 1e-6 * lib.size
    cpms <- x$count / lib.size
  }

  cpms
}

##' @importFrom edgeR rpkm
##' @export
rpkm <- function(x, ...) UseMethod("rpkm")

##' Cacluate RPKM from a FacileDataSet
##'
##' @method rpkm tbl_sqlite
##' @rdname rpkm
##' @importFrom edgeR rpkm
##' @export
##' @param x an expression-like facile result
##' @param lib.size ignored for now, this is fetched from the
##'   \code{FacileDataSet}
##' @param log log the result?
##' @param prior.count prior.count to add to observed counts
##' @param .fds A \code{FacileDataSet} object
##' @return a modified expression-like result with a \code{rpkm} column.
rpkm.tbl_sqlite <- function(x, gene.length=NULL, lib.size=NULL, log=FALSE,
                            prior.count=5, sample.stats=NULL,
                            .fds=fds(x), ...) {
  stopifnot(is.FacileDataSet(.fds))
  if (is.null(gene.length)) {
    gene.length <- gene_info_tbl(.fds) %>% collect
  }
  stopifnot(all(c('feature_id', 'length') %in% colnames(gene.length)))
  cpms <- cpm(x, lib.size=lib.size, log=log, prior.count=prior.count,
              sample.stats=sample.stats, .fds=.fds, ...)
  out <- calc.rpkm(cpms, gene.length, log)
  set_fds(out, .fds)
}

##' @method rpkm tbl_sqlite
##' @importFrom edgeR rpkm
##' @export
##' @rdname rpkm
rpkm.tbl_df <- function(x, gene.length=NULL, lib.size=NULL, log=FALSE,
                        prior.count=5, sample.stats=NULL, .fds=fds(x), ...) {
  assert_expression_result(x)
  stopifnot(is.FacileDataSet(.fds))
  if (is.null(gene.length)) {
    gene.length <- gene_info_tbl(.fds) %>% collect
  }
  stopifnot(all(c('feature_id', 'length') %in% colnames(gene.length)))
  cpms <- cpm(x, lib.size=lib.size, log=log, prior.count=prior.count,
              sample.stats=sample.stats, .fds=.fds, ...)
  out <- calc.rpkm(cpms, gene.length, log=log)
  set_fds(out, .fds)
}

##' @method rpkm tbl_dt
##' @export
rpkm.tbl_dt <- function(x, gene.length=NULL, lib.size=NULL, log=FALSE,
                        prior.count=5, sample.stats=NULL, .fds=fds(x), ...) {
  rpkm.tbl_df(x, gene.length, lib.size, log, prior.count, sample.stats,
              .fds, ...)
}

calc.rpkm <- function(cpms, gene.length, log) {
  assert_expression_result(cpms)
  stopifnot('cpm' %in% colnames(cpms))
  stopifnot(all(c('feature_id', 'length') %in% colnames(gene.length)))
  xref <- match(cpms[['feature_id']], gene.length[['feature_id']])
  kb <- gene.length[['length']][xref] / 1000
  if (log) {
    cpms[['rpkm']] <- cpms[['cpm']] - log2(kb)
  } else {
    cpms[['rpkm']] <- cpms[['cpm']] / kb
  }
  cpms
}

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
##' @param covariates A \code{character} vector specifying the  additional
##'   covariates to append to \code{out$samples}. Must be valid entries in the
##'   \code{sample_covariate::variable} column.
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

##' @export
##' @method as.DGEList matrix
##' @rdname expression-container
as.DGEList.matrix <- function(x, covariates=NULL, feature_ids=NULL,
                              .fds=fds(x), custom_key=NULL, ...) {
  stopifnot(is(x, 'FacileExpression'))
  .fds <- force(.fds)
  stopifnot(is.FacileDataSet(.fds))

  samples <- tibble(
    dataset=sub('_.*$', '', colnames(x)),
    sample_id=sub('^.*?_', '', colnames(x)))

  fids <- rownames(x)
  genes <- gene_info_tbl(.fds) %>%
    filter(feature_id %in% fids) %>%
    collect(n=Inf) %>%
    as.data.frame %>%
    set_rownames(., .$feature_id)

  class(x) <- 'matrix'
  y <- DGEList(x, genes=genes)

  ## Doing the internal filtering seems to be too slow
  ## sample.stats <- fetch_sample_statistics(db, x) %>%
  sample.stats <- fetch_sample_statistics(.fds, samples) %>%
    collect(n=Inf) %>%
    mutate(samid=paste(dataset, sample_id, sep='_')) %>%
    rename(lib.size=libsize, norm.factors=normfactor) %>%
    as.data.frame %>%
    set_rownames(., .$samid)

  y$samples <- cbind(
    transform(y$samples, lib.size=NULL, norm.factors=NULL),
    sample.stats[colnames(y),,drop=FALSE])

  if (!is.null(covariates)) {
    if (is.character(covariates) && length(covariates)) {
      covariates <- fetch_sample_covariates(.fds, samples, covariates,
                                            custom_key)
    }
    assert_sample_covariates(covariates)
    covs <- spread_covariates(covariates, .fds) %>%
      as.data.frame %>%
      set_rownames(., paste(.$dataset, .$sample_id, sep='_')) %>%
      select(-dataset, -sample_id)
    y$samples <- cbind(y$samples, covs[colnames(y),,drop=FALSE])
  }

  if (!is.null(feature_ids) && is.character(feature_ids)) {
    keep <- feature_ids %in% rownames(y)
    if (mean(keep) != 1) {
      warning(sprintf("Only %d / %d feature_ids requested are in dataset",
                      sum(keep), length(keep)))
    }
    y <- y[feature_ids[keep],]
  }
  set_fds(y, .fds)
}

##' @export
##' @method as.DGEList data.frame
##' @rdname expression-container
as.DGEList.data.frame <- function(x, covariates=NULL, feature_ids=NULL,
                                  .fds=fds(x), custom_key=NULL, ...) {
  .fds <- force(.fds)
  stopifnot(is.FacileDataSet(.fds))

  x <- assert_sample_subset(x)
  has.count <- 'count' %in% colnames(x)
  refetch <- is.null(feature_ids) && !has.count

  if (refetch) {
    if (has.count) {
      warning("Ignoring expression in `x` and fetching data for `feature_ids`",
              immediate.=TRUE)
    }
    counts <- fetch_expression(.fds, x, feature_ids=feature_ids, as.matrix=TRUE)
  }  else {
    counts.dt <- assert_expression_result(x) %>%
      collect(n=Inf) %>%
      setDT %>%
      unique(by=c('dataset', 'sample_id', 'feature_id'))
    counts.dt[, samid := paste(dataset, sample_id, sep='_')]
    counts <- local({
      wide <- dcast.data.table(counts.dt, feature_id ~ samid, value.var='count')
      out <- as.matrix(wide[, -1, with=FALSE])
      rownames(out) <- wide[[1L]]
      class(out) <- c('FacileExpression', class(out))
      out
    })
  }

  as.DGEList(counts, covariates=covariates, feature_ids=feature_ids,
             .fds=.fds, custom_key=custom_key, ...)
}

##' @export
##' @rdname expression-container
as.ExpressionSet <- function(x, covariates=NULL, feature_ids=NULL,
                             exprs='counts', .fds=fds(x), custom_key=NULL,
                             ...) {
  .fds <- force(.fds)
  stopifnot(is.FacileDataSet(.fds))
  assert_sample_subset(x)
  if (!require("Biobase")) stop("Biobase required")
  y <- as.DGEList(x, covariates, feature_ids, .fds=.fds, custom_key=custom_key,
                  ...)
  es <- ExpressionSet(y$counts)
  pData(es) <- y$samples
  fData(es) <- y$genes
  set_fds(es, .fds)
}
