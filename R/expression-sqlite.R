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
##'
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
fetch_expression.db <- function(x, samples=NULL, feature_ids=NULL,
                                with_symbols=TRUE, do.collect=FALSE) {
  stopifnot(is.FacileDataSet(x))
  dat <- expression_tbl(x)

  if (!is.null(feature_ids)) {
    assertCharacter(feature_ids)
    feature_ids <- unique(feature_ids)

    # if (length(feature_ids) == 1L) {
    #   dat <- filter(dat, feature_id == feature_ids)
    # } else if (length(feature_ids) > 1L) {
    #   dat <- filter(dat, feature_id %in% feature_ids)
    # }

    if (length(feature_ids) == 1L) {
      genes <- gene_info_tbl(x) %>% filter(feature_id == feature_ids)
    } else if (length(feature_ids) > 1L) {
      genes <- gene_info_tbl(x) %>% filter(feature_id %in% feature_ids)
    }

    if (with_symbols) {
      genes %<>% select(feature_id, symbol)
    } else {
      genes %<>% select(feature_id)
    }

    dat <- inner_join(dat, genes, by='feature_id')
  }

  dat <- filter_samples(dat, samples)

  if (do.collect) {
    dat <- collect(dat, n=Inf)
  }

  class(dat) <- c('FacileExpression', class(dat))
  set_fds(dat, x)
}

##' Append expression values to sample-descriptor
##'
##' @export
##' @param x a samples descriptor
##' @param feature_ids character vector of feature_ids
##' @param with_symbols Do you want gene symbols returned, too?
##' @param .fds A \code{FacileDataSet} object
##' @return a tbl-like result
with_expression.db <- function(samples, feature_ids, with_symbols=TRUE,
                               .fds=fds(samples)) {
  stopifnot(is.FacileDataSet(.fds))
  stopifnot(is.character(feature_ids) && length(feature_ids) > 0)
  samples <- assert_sample_subset(samples)

  .fds %>%
    fetch_expression.db(samples, feature_ids, with_symbols=with_symbols) %>%
    join_samples(samples) %>%
    set_fds(.fds)
}

##' @export
cpmdb <- function(x, ...) UseMethod("cpmdb")

##' Calculated counts per million from a FacileDataSet
##'
##' @method cpmdb tbl_sqlite
##' @rdname cpmdb
##'
##' @export
##' @importFrom edgeR cpm
##'
##' @param x an expression-like facile result
##' @param lib.size ignored for now, this is fetched from the
##'   \code{FacileDataSet}
##' @param log log the result?
##' @param prior.count prior.count to add to observed counts
##' @param .fds A \code{FacileDataSet} object
##' @return a modified expression-like result with a \code{cpm} column.
cpmdb.tbl_sqlite <- function(x, lib.size=NULL, log=FALSE, prior.count=5,
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

##' @method cpmdb tbl_df
##' @rdname cpmdb
##' @export
cpmdb.tbl_df <- function(x, lib.size=NULL, log=FALSE, prior.count=5,
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

  tmp <- inner_join(x, sample.stats, by=c('dataset', 'sample_id'))
  cpms <- calc.cpm(tmp, lib.size=tmp$libsize * tmp$normfactor, log=log,
                   prior.count=prior.count)
  tmp %>%
    mutate(cpm=cpms) %>%
    select_(.dots=c(colnames(x), 'cpm')) %>%
    set_fds(.fds)
}

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

##' @export rpkmdb
rpkmdb <- function(x, ...) UseMethod("cpmdb")

##' Cacluate RPKM from a FacileDataSet
##'
##' @method rpkmdb tbl_sqlite
##' @rdname rpkmdb
##' @importFrom edgeR rpkm
##' @export
##' @param x an expression-like facile result
##' @param lib.size ignored for now, this is fetched from the
##'   \code{FacileDataSet}
##' @param log log the result?
##' @param prior.count prior.count to add to observed counts
##' @param .fds A \code{FacileDataSet} object
##' @return a modified expression-like result with a \code{rpkm} column.
rpkmdb.tbl_sqlite <- function(x, gene.length=NULL, lib.size=NULL, log=FALSE,
                            prior.count=5, .fds=fds(x), ...) {
  stopifnot(is.FacileDataSet(.fds))
  if (is.null(gene.length)) {
    gene.length <- gene_info_tbl(.fds) %>% collect
  }
  stopifnot(all(c('feature_id', 'length') %in% colnames(gene.length)))
  cpms <- cpm(x, lib.size=lib.size, log=log, prior.count=prior.count,
              .fds=.fds, ...)
  out <- calc.rpkm(cpms, gene.length, log)
  set_fds(out, .fds)
}

##' @method rpkmdb tbl_sqlite
##' @importFrom edgeR rpkm
##' @export
##' @rdname rpkmdb
rpkmdb.tbl_df <- function(x, gene.length=NULL, lib.size=NULL, log=FALSE,
                        prior.count=5, .fds=fds(x), ...) {
  assert_expression_result(x)
  stopifnot(is.FacileDataSet(.fds))
  if (is.null(gene.length)) {
    gene.length <- gene_info_tbl(.fds) %>% collect
  }
  stopifnot(all(c('feature_id', 'length') %in% colnames(gene.length)))
  cpms <- cpm(x, lib.size=lib.size, log=log, prior.count=prior.count,
              .fds=.fds, ...)
  out <- calc.rpkm(cpms, gene.length, log=log)
  set_fds(out, .fds)
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

##' Converts a result from `fetch_expression` into a DGEList
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
##' @param .fds The \code{FacileDataSet} that \code{x} was retrieved from.
##' @return a \code{\link[edgeR]{DGEList}}
as.DGEList.db <- function(x, covariates=NULL, .fds=fds(x), ...) {
  .fds <- force(.fds)
  stopifnot(is.FacileDataSet(.fds))

  counts.df <- assert_expression_result(x) %>%
    collect(n=Inf) %>%
    distinct(dataset, sample_id, feature_id, .keep_all=TRUE) %>%
    mutate(samid=paste(dataset, sample_id, sep='_'))
  counts <- acast(counts.df, feature_id ~ samid, value.var='count')

  samples <- counts.df %>%
    select(dataset, sample_id) %>%
    distinct(.keep_all=TRUE)

  fids <- rownames(counts)
  genes <- gene_info_tbl(.fds) %>%
    filter(feature_id %in% fids) %>%
    collect(n=Inf) %>%
    as.data.frame %>%
    set_rownames(., .$feature_id)

  y <- DGEList(counts, genes=genes)

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

  if (is.character(covariates) && length(covariates)) {
    covs <- with_sample_covariates(samples, covariates, .fds) %>%
      as.data.frame %>%
      set_rownames(., paste(.$dataset, .$sample_id, sep='_')) %>%
      select(-dataset, -sample_id)
    y$samples <- cbind(y$samples, covs[colnames(y),,drop=FALSE])
  }

  set_fds(y, .fds)
}

##' Create an ExpressionSet from `fetch_expression`.
##'
##' @rdname expression-container
##' @export
##' @param x a facile expression-like result
##' @param covariates A \code{character} vector specifying the  additional
##'   covariates to append to \code{out$samples}. Must be valid entries in the
##'   \code{sample_covariate::variable} column.
##' @param assay Which column to put in \code{"exprs"}
##' @param .fds The \code{FacileDataSet} that \code{x} was retrieved from.
##' @return a \code{\link[Biobase]{ExpressionSet}}
as.ExpressionSet <- function(x, covariates=NULL, exprs='counts', .fds=fds(x),
                             ...) {
  if (!require("Biobase")) {
    stop("Biobase required")
  }
  y <- as.DGEList(x, covariates, db, .fds=.fds, ...)
  es <- ExpressionSet(y$counts)
  pData(es) <- y$samples
  fData(es) <- y$genes
  set_fds(es, .fds)
}
