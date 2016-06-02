##' Helper function creates dplyr query to get expression data of interest.
##'
##' Note that if \code{db} is not provided, the result will be \code{collect}ed,
##' because the database connection will be closed when the function exits.
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
    dat <- collect(dat)
    db <- NULL
  }

  class(dat) <- c('FacileExpression', class(dat))
  attr(dat, 'db') <- db
  dat
}

##' Append expression values to sample-descriptor
##'
##' @export
##' @param x a samples descriptor
##' @param feature_ids character vector of feature_ids
##' @param a \code{FacileDb} object
##' @return a tbl-like result
with_expression <- function(samples, feature_ids, db=attr(samples, 'db')) {
  stopifnot(is.FacileDb(db))
  stopifnot(is.character(feature_ids) && length(feature_ids) > 0)
  samples <- assert_sample_subset(samples)

  out <- fetch_expression(db, samples, feature_ids) %>%
    join_samples(samples)
  attr(out, 'db') <- db
  out
}

##' @method cpm tbl_sqlite
##' @export
cpm.tbl_sqlite <- function(x, lib.size=NULL, log=FALSE, prior.count=5,
                           db=attr(x, 'db'), ...) {
  stopifnot(is.FacileDb(db))
  ## Let's leverage the fact that we're already working in the database
  sample.stats <- fetch_sample_statistics(db, x)
  out <- cpm(collect(x), lib.size=lib.size, log=log, prior.count=prior.count,
             sample.stats=sample.stats, db=db, ...)
  attr(out, 'db') <- db
  out
}

##' @method cpm tbl_df
##' @importFrom edgeR cpm
##' @export
cpm.tbl_df <- function(x, lib.size=NULL, log=FALSE, prior.count=5,
                       sample.stats=NULL, db=attr(x, 'db'), ...) {
  assert_expression_result(x)
  stopifnot(is.FacileDb(db))
  if (is.null(sample.stats)) {
    sample.stats <- fetch_sample_statistics(x, db=db)
  }
  sample.stats <- sample.stats %>%
    assert_sample_statistics %>%
    collect

  ## cpm.DGList functionality unrolled over a data.frame
  mult <- mutate(sample.stats, lib.size=libsize * normfactor,
                 pr.count=lib.size / mean(lib.size) * prior.count)

  xref <- match(paste0(x$dataset, x$sample_id),
                paste0(mult$dataset, mult$sample_id))

  if (log) {
    mult$lib.size <- 1e-6 * (mult$lib.size + 2*mult$pr.count)
    x$cpm <- log2((x$count + mult$pr.count[xref]) / mult$lib.size[xref])
  } else {
    mult$lib.size <- 1e-6 * mult$lib.size
    x$cpm <- x$count / mult$lib.size[xref]
  }

  x
}


##' Converts a result from `fetch_expression` into a DGEList
##'
##' The genes and samples that populate the \code{DGEList} are specified by
##' \code{x}, and the caller can request addition sample information to be
##' appended to \code{out$samples} via specification through the
##' \code{covariates} argument.
##'
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
as.DGEList <- function(x, covariates=NULL, db=attr(x, 'db'),
                       cov.def=db[['cov.def']], ...) {
  if (FALSE) {
    covariates <- c('IC', 'TC', 'BCOR')
  }
  stopifnot(is.FacileDb(db))
  counts <- assert_expression_result(x) %>%
    collect %>%
    mutate(samid=paste(dataset, sample_id, sep='_')) %>%
    mcast(feature_id ~ samid, value.var='count')

  fids <- rownames(counts)
  genes <- gene_info_tbl(db) %>%
    filter(feature_id %in% fids) %>%
    collect %>%
    as.data.frame %>%
    set_rownames(., .$feature_id)

  y <- edgeR::DGEList(counts, genes=genes)

  ## Doing the internal filtering seems to be too slow
  ## sample.stats <- fetch_sample_statistics(db, x) %>%
  sample.stats <- fetch_sample_statistics(db) %>%
    collect %>%
    mutate(samid=paste(dataset, sample_id, sep='_')) %>%
    rename(lib.size=libsize, norm.factors=normfactor)

  xref <- match(colnames(y), sample.stats$samid)
  # y$samples$lib.size <- sample.stats$libsize[xref]
  # y$samples$norm.factors <- sample.stats$normfactor[xref]
  # y$samples <- cbind(y$samples, sample.stats[xref, c('dataset', 'sample_id')])

  y$samples <- cbind(
    transform(y$samples, lib.size=NULL, norm.factors=NULL),
    sample.stats[xref, c('lib.size', 'norm.factors', 'dataset', 'sample_id')])

  if (is.character(covariates)) {
    covs <- fetch_sample_covariates(db, y$samples, covariates) %>%
      spread_covariates(cov.def) %>%
      select(-dataset, -sample_id) %>%
      extract(colnames(y),,drop=FALSE) %>%
      set_rownames(colnames(y)) ## to fix NA, NA.1, rownames if missing samples
    y$samples <- cbind(y$samples, covs)
  }

  y
}
