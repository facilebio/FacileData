##' Utility functions to get row and column indices of rnaseq hdf5 files.
##'
##' @export
##' @rdname hdf5expression
hdf5_sample_indices <- function(x, assay_name='rnaseq', samples=NULL) {
  stopifnot(is.FacileDataSet(x))
  assert_string(assay_name)
  stopifnot(assay_name %in% assay_names(x))

  asi <- assay_sample_info_tbl(x) %>%
    filter(assay == assay_name) %>%
    select(assay, dataset, sample_id, hdf5_index) %>%
    collect(n=Inf)

  if (is.null(samples)) {
    samples <- asi
  } else {
    samples <- samples %>%
      assert_sample_subset %>%
      distinct(dataset, sample_id) %>%
      collect(n=Inf) %>%
      mutate(assay=assay_name) %>%
      left_join(asi, by=c('assay', 'dataset', 'sample_id'))
  }

  samples
}

##' @export
##' @rdname hdf5expression
hdf5_gene_indices <- function(x, assay_name='rnaseq', feature_ids=NULL) {
  stopifnot(is.FacileDataSet(x))
  if (is.null(feature_ids)) {
    feature_ids <- character()
  }
  assert_character(feature_ids)
  feature_ids <- unique(feature_ids)

  fi <- feature_info_tbl(x) %>%
    filter(feature_type == 'entrez')
  if (length(feature_ids) == 1) {
    fi <- filter(fi, feature_id == feature_ids)
  } else if (length(feature_ids) > 1) {
    fi <- filter(fi, feature_id %in% feature_ids)
  }
  fi <- collect(fi, n=Inf)

  genes <- assay_feature_info_tbl(x) %>%
    filter(assay == assay_name) %>%
    collect(n=Inf) %>%
    inner_join(fi, by='feature_id') %>%
    rename(symbol=name)

  genes
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
##' @param samples a sample descriptor to specify which samples to return data
##'   from.
##' @param feature_ids character vector of ids, if \code{NULL} returns all
##'   features
##' @param with_symbols adds a symbol column to the outgoing result
##' @param as.matrix by default, the data is returned in a long-form tbl-like
##'   result. If set to \code{TRUE}, the data is returned as a matrix.
##' @return A lazy \code{\link[dplyr]{tbl}} object with the expression
##'   data to be \code{\link[dplyr]{collect}}ed when \code{db} is provided,
##'   othwerise a \code{tbl_df} of the results.
fetch_expression <- function(x, samples=NULL, feature_ids=NULL,
                             with_symbols=!as.matrix, as.matrix=FALSE) {
  ## In the 'unhinged'/multiassay we force assay='rnaseq'
  stopifnot(is.FacileDataSet(x))

  hdf.sample.idxs <- hdf5_sample_indices(x, 'rnaseq', samples)

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
  ask.some <- is.character(feature_ids)
  fetch.all <- !ask.some || length(feature_ids) > 700

  gene.info <- hdf5_gene_indices(x, feature_ids=feature_ids)
  gene.info <- arrange(gene.info, hdf5_index)

  ridx <- if (fetch.all) NULL else gene.info$hdf5_index

  if (isTRUE(as.matrix) && isTRUE(with_symbols)) {
    warning("with_symbols ignored when as.matrix=TRUE")
  }

  counts <- hdf.sample.idxs %>%
    group_by(dataset) %>%
    do(res={
      ds <- .$dataset[1L]
      hd5.name <- paste0('assay/rnaseq/', ds)
      cnts <- h5read(x$hdf5.fn, hd5.name, list(ridx, .$hdf5_index))
      stopifnot(nrow(cnts) == nrow(gene.info))
      if (ask.some && fetch.all) {
        cnts <- cnts[gene.info$hdf5_index,,drop=FALSE]
      }
      if (isTRUE(as.matrix)) {
        dimnames(cnts) <- list(gene.info$feature_id,
                               paste(ds, .$sample_id, sep='_'))
      } else {
        dimnames(cnts) <- list(gene.info$feature_id, .$sample_id)
        cnts <- as.data.table(cnts, keep.rownames=TRUE)
        cnts <- melt.data.table(cnts, id.vars='rn')
        cnts[, dataset := ds]
        cnts[, variable := as.character(variable)]
        setnames(cnts, c('feature_id', 'sample_id', 'count', 'dataset'))
        setcolorder(cnts, c('dataset', 'sample_id', 'feature_id', 'count'))
        if (isTRUE(with_symbols)) {
          sxref <- match(cnts$feature_id, gene.info$feature_id)
          cnts[, symbol := gene.info$symbol[sxref]]
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

