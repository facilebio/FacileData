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
##' @param db The db connection to use. Creates and destroys a new one
##'   if not passed in.
##' @param indication The indication to get data from
##' @param subtype The subtype from the desired \code{indication}
##' @param sample_ids The sample_id's (barcodes) to fetch. If specified
##'   this trumps whatever is specified in \code{indication} and
##'   \code{subtype}
##' @param entrez_ids Restricts the fetch to only these genes
##' @return A lazy \code{\link[dplyr]{tbl}} object with the expression
##'   data to be \code{\link[dplyr]{collect}}ed when \code{db} is provided,
##'   othwerise a \code{tbl_df} of the results.
fetch_expression <- function(db, samples=NULL, feature_ids=NULL,
                             do.collect=FALSE) {
  stopifnot(is.FacileDb(db))
  dat <- expression_tbl(db)

  if (!is.null(feature_ids)) {
    assertCharacter(feature_ids)
    feature_ids <- unique(feature_ids)
    if (length(feature_ids) == 1L) {
      dat <- filter(dat, feature_id == feature_ids)
    } else if (length(feature_ids) > 1L) {
      dat <- filter(dat, feature_id %in% feature_ids)
    }
  }

  dat <- filter_samples(dat, samples)

  if (do.collect) {
    dat <- collect(dat, n=Inf)
    db <- NULL
  }

  class(dat) <- c('FacileExpression', class(dat))
  set_fdb(dat, db)
}

##' Append expression values to sample-descriptor
##'
##' @export
##' @param x a samples descriptor
##' @param feature_ids character vector of feature_ids
##' @param a \code{FacileDb} object
##' @return a tbl-like result
with_expression <- function(samples, feature_ids, db=fdb(samples)) {
  stopifnot(is.FacileDb(db))
  stopifnot(is.character(feature_ids) && length(feature_ids) > 0)
  samples <- assert_sample_subset(samples)

  out <- fetch_expression(db, samples, feature_ids) %>%
    join_samples(samples)
  set_fdb(out, db)
}

##' @method cpm tbl_sqlite
##' @importFrom edgeR cpm
##' @export
cpm.tbl_sqlite <- function(x, lib.size=NULL, log=FALSE, prior.count=5,
                           db=fdb(x), ...) {
  stopifnot(is.FacileDb(db))
  if (!is.null(lib.size)) {
    warning("not supporting custom lib.size yet")
  }
  ## Let's leverage the fact that we're already working in the database
  sample.stats <- fetch_sample_statistics(db)

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
               prior.count=prior.count, sample.stats=sample.stats, db=db, ...)

  }
  set_fdb(out, db)
}

##' @method cpm tbl_df
##' @export
cpm.tbl_df <- function(x, lib.size=NULL, log=FALSE, prior.count=5,
                       sample.stats=NULL, db=fdb(x), ...) {
  assert_expression_result(x)
  stopifnot(is.FacileDb(db))
  if (!is.null(lib.size)) {
    warning("not supporting custom lib.size yet")
  }

  if (is.null(sample.stats)) {
    samples <- distinct(x, dataset, sample_id)
    sample.stats <- fetch_sample_statistics(db, samples)
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
    set_fdb(db)
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

##' @method rpkm tbl_sqlite
##' @importFrom edgeR rpkm
##' @export
rpkm.tbl_sqlite <- function(x, gene.length=NULL, lib.size=NULL, log=FALSE,
                            prior.count=5, db=fdb(x), ...) {
  stopifnot(is.FacileDb(db))
  if (is.null(gene.length)) {
    gene.length <- gene_info_tbl(db) %>% collect
  }
  stopifnot(all(c('feature_id', 'length') %in% colnames(gene.length)))
  out <- cpm(x, lib.size=lib.size, log=log, prior.count=prior.count, db=db, ...)
  calc.rpkm(out, gene.length, log=log)
}

##' @method rpkm tbl_sqlite
##' @importFrom edgeR rpkm
##' @export
rpkm.tbl_df <- function(x, gene.length=NULL, lib.size=NULL, log=FALSE,
                        prior.count=5, db=fdb(x), ...) {
  assert_expression_result(x)
  stopifnot(is.FacileDb(db))
  if (is.null(gene.length)) {
    gene.length <- gene_info_tbl(db) %>% collect
  }
  stopifnot(all(c('feature_id', 'length') %in% colnames(gene.length)))
  out <- cpm(x, lib.size=lib.size, log=log, prior.count=prior.count, db=db, ...)
  calc.rpkm(out, gene.length, log=log)
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
##' @param x \code{tbl_sql} result from a call to \code{\link{fetch_expression}}
##' @param covariates A \code{character} vector specifying the  additional
##'   covariates to append to \code{out$samples}. Must be valid entries in the
##'   \code{sample_covariate::variable} column.
##' @param db The \code{FacilDb} object. This is extracted from \code{x} if
##'   we are able.
##' @param cov.def the path to the yaml file that defines what each type of
##'   variable is. This is also set in and extracted from \code{db}.
##' @return a \code{\link[edgeR]{DGEList}}
as.DGEList <- function(x, covariates=NULL, db=fdb(x),
                       cov.def=db[['cov.def']], ...) {
  if (FALSE) {
    covariates <- c('stage', 'sex')
  }
  stopifnot(is.FacileDb(db))

  counts.df <- assert_expression_result(x) %>%
    collect(n=Inf) %>%
    distinct(dataset, sample_id, feature_id, .keep_all=TRUE) %>%
    mutate(samid=paste(dataset, sample_id, sep='_'))
  counts <- mcast(counts.df, feature_id ~ samid, value.var='count')

  samples <- counts.df %>%
    select(dataset, sample_id) %>%
    distinct(.keep_all=TRUE)

  fids <- rownames(counts)
  genes <- gene_info_tbl(db) %>%
    filter(feature_id %in% fids) %>%
    collect(n=Inf) %>%
    as.data.frame %>%
    set_rownames(., .$feature_id)

  y <- edgeR::DGEList(counts, genes=genes)

  ## Doing the internal filtering seems to be too slow
  ## sample.stats <- fetch_sample_statistics(db, x) %>%
  sample.stats <- fetch_sample_statistics(db, samples) %>%
    collect(n=Inf) %>%
    mutate(samid=paste(dataset, sample_id, sep='_')) %>%
    rename(lib.size=libsize, norm.factors=normfactor) %>%
    as.data.frame %>%
    set_rownames(., .$samid)

  y$samples <- cbind(
    transform(y$samples, lib.size=NULL, norm.factors=NULL),
    sample.stats[colnames(y),,drop=FALSE])

  if (is.character(covariates) && length(covariates)) {
    covs <- with_sample_covariates(samples, covariates, db) %>%
      as.data.frame %>%
      set_rownames(., paste(.$dataset, .$sample_id, sep='_')) %>%
      select(-dataset, -sample_id)
    y$samples <- cbind(y$samples, covs[colnames(y),,drop=FALSE])
  }

  y
}

##' Create an ExpressionSet from `fetch_expression`.
##'
##' @rdname expression-container
##' @export
##' @param x \code{tbl_sql} result from a call to \code{\link{fetch_expression}}
##' @param covariates A \code{character} vector specifying the  additional
##'   covariates to append to \code{out$samples}. Must be valid entries in the
##'   \code{sample_covariate::variable} column.
##' @param assay Which column to put in \code{"exprs"}
##' @param db The \code{FacilDb} object. This is extracted from \code{x} if
##'   we are able.
##' @param cov.def the path to the yaml file that defines what each type of
##'   variable is. This is also set in and extracted from \code{db}.
##' @return a \code{\link[Biobase]{ExpressionSet}}
as.ExpressionSet <- function(x, covariates=NULL, exprs='counts', db=fdb(x),
                             cov.def=db[['cov.def']], ...) {
  if (!require("Biobase")) {
    stop("Biobase required")
  }
  y <- as.DGEList(x, covariates, db, cov.df, ...)
  es <- ExpressionSet(y$counts)
  pData(es) <- y$samples
  fData(es) <- y$genes
  es
}
