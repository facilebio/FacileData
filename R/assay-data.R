## TODO: Investigate why fetching a one gene DGEList across TCGA gives
##       "library size of zero detected" warning, ie.
##       y <- as.DGEList(sample_stats_tbl(fds), feature_id="919", covariates=NULL)
##' Utility functions to get row and column indices of rnaseq hdf5 files.
##'
##' @export
assay_sample_info <- function(x, assay_name, samples=NULL) {
  stopifnot(is.FacileDataSet(x))
  if (!is.null(samples)) {
    samples <- assert_sample_subset(samples) %>%
      collect(n=Inf) %>%
      distinct(dataset, sample_id)
  }
  feature.type <- assay_feature_type(x, assay_name) ## validate assay_name

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


##' Returns the feature_type for a given assay
##'
##' The elements of the rows for a given assay all correspond to a particular
##' feature space (ie. feature_type='entrez')
##'
##' @export
##' @param x \code{FacileDataSet}
##' @param assay_name the name of the assay
assay_feature_type <- function(x, assay_name) {
  stopifnot(is.FacileDataSet(x))
  assert_string(assay_name)
  assert_choice(assay_name, assay_types(x))
  assay_info_tbl(x) %>%
    filter(assay == assay_name) %>%
    collect %$%
    feature_type
}

##' Materializes a table with all feature information for a given assay.
##'
##' This fetches the hdf5_index for the assays as well
##' @export
##' @inheritParams assay_feature_type
##' @param feature_ids a character vector of feature_ids
##' @return a \code{tbl_sqlite} result with the feature information for the
##'   features in a specified assay. This
assay_feature_info <- function(x, assay_name, feature_ids=NULL) {
  ftype <- assay_feature_type(x, assay_name)
  if (!is.null(feature_ids)) {
    assert_character(feature_ids)
    if (length(feature_ids) == 0) {
      feature_ids <- NULL
    } else {
      feature_ids <- tibble(assay=assay_name, feature_id=feature_ids)
    }
  }

  if (is.null(feature_ids)) {
    afinfo <- assay_feature_info_tbl(x) %>% filter(assay == assay_name)
  } else {
    afinfo <- assay_feature_info_tbl(x) %>%
      inner_join(feature_ids, by=c('assay', 'feature_id'), copy=TRUE,
                 auto_index=TRUE)
  }

  assay.info <- select(assay_info_tbl(x), assay, assay_type, feature_type)
  out <- afinfo %>%
    inner_join(assay.info, by='assay') %>%
    inner_join(feature_info_tbl(x), by=c('feature_type', 'feature_id'))
  set_fds(out, x)
}


##' Fetch data from single assay of choice
##'
##' @export
##' @importFrom rhdf5 h5read
##' @inheritParams assay_feature_info
##' @param x A \code{FacileDataSet} object.
##' @param samples a sample descriptor to specify which samples to return data
##'   from.
##' @param normalized return normalize or raw data values, defaults to
##'   \code{raw}
##' @param as.matrix by default, the data is returned in a long-form tbl-like
##'   result. If set to \code{TRUE}, the data is returned as a matrix.
##' @param ... parameters to pass to normalization methods
##' @return A lazy \code{\link[dplyr]{tbl}} object with the expression
##'   data to be \code{\link[dplyr]{collect}}ed when \code{db} is provided,
##'   othwerise a \code{tbl_df} of the results.
fetch_assay_data <- function(x, assay_name, feature_ids=NULL, samples=NULL,
                             normalized=TRUE, as.matrix=FALSE, ...,
                             .subset.threshold=700) {
  stopifnot(is.FacileDataSet(x))
  assert_flag(as.matrix)
  finfo <- assay_feature_info(x, assay_name, feature_ids=feature_ids) %>%
    collect(n=Inf) %>%
    arrange(hdf5_index)
  atype <- finfo$assay_type[1L]
  ftype <- finfo$feature_type[1L]
  sinfo <- assay_sample_info(x, assay_name, samples) %>%
    mutate(samid=paste(dataset, sample_id, sep="_"))

  ## DEBUG: Tune chunk size?
  ## As the number of genes you are fetching increases, only subsetting
  ## out a few of them intead of first loading the whole matrix gives you little
  ## value, for instance, over all BRCA tumors (994) these are some timings for
  ## different numbers of genes:
  ##
  ##     Ngenes                   Time (seconds)
  ##     10                       0.5s
  ##     100                      0.8s
  ##     250                      1.2s
  ##     500                      2.6s
  ##     750                      6s seconds
  ##     3000                     112 seconds!
  ##     unpspecified (all 26.5k) 7 seconds!
  ##
  ## I'm using `ridx` as a hack downstream to run around the issue of slowing
  ## down the data by trying to subset many rows via hdf5 instead of loading
  ## then subsetting after (this is so weird)
  ##
  ## TODO: setup unit tests to ensure that ridx subsetting and remapping back
  ## to original genes works
  fetch.some <- nrow(finfo) < .subset.threshold
  ridx <- if (fetch.some) finfo$hdf5_index else NULL

  dat <- sinfo %>%
    group_by(dataset) %>%
    do(res={
      ds <- .$dataset[1L]
      hd5.name <- paste('assay', assay_name, ds, sep='/')
      vals <- h5read(hdf5fn(x), hd5.name, list(ridx, .$hdf5_index))
      if (is.null(ridx)) {
        vals <- vals[finfo$hdf5_index,,drop=FALSE]
      }
      dimnames(vals) <- list(finfo$feature_id, .$samid)
      if (normalized) {
        vals <- normalize.assay.matrix(vals, finfo, ., ...)
      }
      if (!as.matrix) {
        vals <- as.data.table(vals, keep.rownames=TRUE)
        vals <- melt.data.table(vals, id.vars='rn', variable.factor=FALSE,
                                variable.name='sample_id')
        setnames(vals, 1L, 'feature_id')
        vals[, dataset := ds]
        vals[, assay := assay_name]
        vals[, assay_type := atype]
        vals[, feature_type := ftype]
        xref <- match(vals$feature_id, finfo$feature_id)
        vals[, feature_name := finfo$name[xref]]
        corder <- c('dataset', 'sample_id', 'assay', 'assay_type',
                    'feature_type', 'feature_id', 'feature_name', 'value')
        setcolorder(vals, corder)
        vals
      }
      vals
    }) %>%
    ungroup

  if (as.matrix) {
    out <- do.call(cbind, dat$res)
  } else {
    out <- as.tbl(setDF(rbindlist(dat$res)))
  }

  class(out) <- c('FacileExpression', class(out))
  set_fds(out, x)
}

## helper functino to fetch_assay_data
normalize.assay.matrix <- function(vals, feature.info, sample.info, ...) {
  stopifnot(
    nrow(vals) == nrow(feature.info),
    all(rownames(vals) == feature.info$feature_id),
    ncol(vals) == nrow(sample.info),
    all(colnames(vals) == sample.info$samid),
    is.character(feature.info$assay_type),
    length(unique(feature.info$assay_type)) == 1L,
    is.numeric(sample.info$libsize), is.numeric(sample.info$normfactor))
  atype <- feature.info$assay_type[1L]
  libsize <- sample.info$libsize * sample.info$normfactor
  if (atype == 'rnaseq') {
    out <- edgeR::cpm(vals, libsize, log=TRUE, prior.count=5)
  } else {
    stop("implement normalization for other assay types")
  }
  out
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

##' Takes a result from fetch_expression and spreads out genes acorss columns
##'
##' This is a convenience function, and will try to guess what you mean if you
##' don't explicitly specify which columns to spread and what to call them.
##' With that mind set, if we find a cpm and symbol column, we will use them
##' because those are the thing you will likely want to use for exploratory
##' data analysis if they're in the incoming dataset. If those columns aren't
##' found, then we'll pick the feature_id and count column.
##'
##' @export
##' @param x facile expression result from \code{fetch_expression}
##' @param key the column from the long-form \code{fetch_expression} table
##'   to put in the columns of the outgoing data.frame that the values are
##'   "spread into"
##' @param value the value column to spread into the \code{key} columns
##' @param .fds the \code{FacileDataSet}
##' @return a more stout \code{x} with the expression values spread across
##'   columns.
spread_expression <- function(x, key=c('symbol', 'feature_id'),
                              value=c('cpm', 'count'), .fds=fds(x)) {
  force(.fds)
  if (missing(key)) {
    key <- if ('symbol' %in% colnames(x)) 'symbol' else 'feature_id'
  }
  if (missing(value)) {
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
                           feature_ids=NULL, as.matrix=FALSE, .fds=fds(x),
                           ...) {
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
               prior.count=prior.count, feature_ids=feature_ids,
               sample.stats=sample.stats, as.matrix=as.matrix, .fds=.fds, ...)

  }
  set_fds(out, .fds)
}

##' @method cpm tbl_df
##' @rdname cpm
##' @export
cpm.tbl_df <- function(x, lib.size=NULL, log=FALSE, prior.count=5,
                       feature_ids=NULL, sample.stats=NULL,
                       as.matrix=as.matrix, .fds=fds(x), ...) {
  stopifnot(is.FacileDataSet(.fds))
  ## TODO: Alter this function to accept `feature_ids` and `as.matrix`
  ## parameter so they can more fluently work in a data piping chain
  # if (!is_expression_result(x)) {
  #   assert_sample_subset(x)
  #   x <- fetch_expression(.fds, x, feature_ids=feature_ids, as.matrix=TRUE)
  # }
  assert_expression_result(x)
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
    gene.length <- gene_info_tbl(.fds) %>% collect(n=Inf)
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
    gene.length <- gene_info_tbl(.fds) %>% collect(n=Inf)
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
                              .fds=fds(x), custom_key=Sys.getenv("USER"), ...) {
  stopifnot(is(x, 'FacileExpression'))
  requireNamespace("edgeR")
  .fds <- force(.fds)
  stopifnot(is.FacileDataSet(.fds))

  ## Construct sample table from colnames of the matrix, and make sure this is
  ## legit
  samples <- tibble(
    dataset=sub('_.*$', '', colnames(x)),
    sample_id=sub('^.*?_', '', colnames(x)))
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
    mutate(samid=paste(dataset, sample_id, sep='_')) %>%
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
      set_rownames(., paste(.$dataset, .$sample_id, sep='_')) %>%
      select(-dataset, -sample_id)
    y$samples <- cbind(y$samples, covs[colnames(y),,drop=FALSE])
  }

  set_fds(y, .fds)
}

##' @export
##' @method as.DGEList data.frame
##' @rdname expression-container
as.DGEList.data.frame <- function(x, covariates=TRUE, feature_ids=NULL,
                                  .fds=fds(x), custom_key=Sys.getenv("USER"),
                                  ...) {
  .fds <- force(.fds)
  stopifnot(is.FacileDataSet(.fds))

  x <- assert_sample_subset(x)

  has.count <- 'count' %in% colnames(x)
  fetch.counts <- !has.count

  ## Do we want to fetch counts from the FacileDataSet?
  if (has.count) {
    if (is.character(feature_ids) && !setequal(feature_ids, x$feature_ids)) {
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
    counts <- fetch_expression(.fds, x, feature_ids=feature_ids, as.matrix=TRUE)
  } else {
    counts.dt <- assert_expression_result(x) %>%
      collect(n=Inf) %>%
      setDT %>%
      unique(by=c('dataset', 'sample_id', 'feature_id'))
    counts.dt[, samid := paste(dataset, sample_id, sep='_')]
    counts <- local({
      wide <- dcast.data.table(counts.dt, feature_id ~ samid, value.var='count')
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
##' @method as.DGEList tbl_sqlite
##' @rdname expression-container
as.DGEList.tbl_sqlite <- function(x, covariates=TRUE, feature_ids=NULL,
                                  .fds=fds(x), custom_key=Sys.getenv("USER"),
                                  ...) {
  x <- collect(x, n=Inf) %>% set_fds(.fds)
  as.DGEList(x, covariates, feature_ids, .fds=.fds, custom_key=custom_key, ...)
}

##' @export
##' @rdname expression-container
as.ExpressionSet <- function(x, covariates=TRUE, feature_ids=NULL,
                             exprs='counts', .fds=fds(x),
                             custom_key=Sys.getenv("USER"), ...) {
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

