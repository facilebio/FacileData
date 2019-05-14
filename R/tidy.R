#' @section DGEList
#' You may have stored extra versions of the count/expression matrix in `x`
#' as a different element of the DGEList. In this case, you can provide a
#' character vector of other expression-like matrices you want to include in
#' the tidying.
#'
#' The following code will include a `"corrected"` column in the tidied output,
#' in addition to the standard `"cpm"` and `"count"` columns.
#'
#' ```
#' y$corrected <- removeBatchEffect(cpm(y), ...)
#' tidy(y, assay_name = "corrected")
#' ```
#'
#' @noRd
#' @export
#' @importFrom edgeR cpm
#' @importFrom reshape2 melt
#' @importFrom broom tidy
#'
#' @method tidy DGEList
tidy.DGEList <- function(x, normalized.lib.sizes = TRUE, prior.count = 3,
                         genes_columns = NULL, samples_columns = NULL,
                         assay_name = NULL, ...) {
  mats <- list(
    cpm = cpm(x, normalized.lib.sizes = normalized.lib.sizes,
              log = TRUE, prior.count = prior.count),
    count = x$counts)

  if (is.character(assay_name)) {
    assay_name <- setdiff(assay_name, "counts")
    for (name in assay_name) {
      m <- x[[name]]
      if (!test_matrix(m, "numeric", nrows = nrow(x), ncols = ncol(x))) {
        stop(sprintf("Cannot extract `%s` matrix from DGEList", name))
      }
      mats[[name]] <- m
    }
  }

  .tidy.core(mats, genes = x$genes, samples = x$samples,
             genes_columns = genes_columns,
             samples_columns = samples_columns, ...)
}

#' @section SummarizedExperiment
#' List as many assay matrices you want to extract from `x` as you like.
#' No further processing is done for these data, so if you want to log2
#' transform the data, you will have to do that first and store it back int
#' `x` as another assay matrix.
#'
#' @noRd
#' @export
#' @importFrom edgeR DGEList
#' @method tidy SummarizedExperiment
tidy.SummarizedExperiment <- function(x, assay_name = NULL,
                                      genes_columns = NULL,
                                      samples_columns = NULL, ...) {
  ns <- tryCatch(loadNamespace("SummarizedExperiment", ...), error = function(e) e)
  if (inherits(ns, "error")) {
    stop("SummarizedExperiment package is required")
  }
  anames <- ns$assayNames(x)
  if (is.null(assay_name)) {
    assay_name <- anames[1L]
  }
  assert_character(assay_name)
  assert_subset(assay_name, anames)

  mats <- lapply(assay_name, function(name) as.matrix(ns$assay(x, name)))
  .tidy.core(mats, genes = ginfo, samples = sinfo,
             genes_columns = genes_columns,
             samples_columns = samples_columns, ...)
}

#' @noRd
#' @export
#' @method tidy EList
tidy.EList <- function(x, assay_name = NULL, ...)  {
  mats <- list(cpm = x$E)
  if (is.matrix(x$weights)) {
    mats$weight <- x$weights
    rownames(mats$weight) <- rownames(x)
    colnames(mats$weight) <- colnames(x)
  } else {
    names(mats)[1L] <- "value"
  }

  if (is.character(assay_name)) {
    assay_name <- setdiff(assay_name, c("E", "weights", "value"))
    for (name in assay_name) {
      m <- x[[name]]
      if (!test_matrix(m, "numeric", nrows = nrow(x), ncols = ncol(x))) {
        stop(sprintf("Cannot extract `%s` matrix from DGEList", name))
      }
      mats[[name]] <- m
    }
  }

  .tidy.core(mats, genes = x$genes, samples = x$targets)
}

#' This is a convenience function for users to whip together a data.frame of
#' assay data annotated with row and column annotations.
#' @noRd
#' @export
tidy.matrix <- function(x, row_covariates, col_covariates, ...) {
  assert_matrix(x, row.names = "unique", col.names = "unique")
  assert_multi_class(row_covariates, c("data.frame", "tbl"))
  assert_multi_class(col_covariates, c("data.frame", "tbl"))
  row_covariates <- as.data.frame(collect(row_covariates, n = Inf))
  col_covariates <- as.data.frame(collect(col_covariates, n = Inf))

  stopifnot(
    nrow(row_covariates) == nrow(x),
    nrow(col_covariates) == ncol(x))

  if (!setequal(rownames(x), rownames(row_covariates))) {
    stop("rownames(x) does not match rownames(row_covariates)")
  }
  row_covariates <- row_covariates[rownames(x),,drop=FALSE]

  if (!setequal(colnames(x), rownames(col_covariates))) {
    stop("colnames(x) does not match rownames(col_covariates)")
  }
  col_covariates <- col_covariates[colnames(x),,drop=FALSE]

  .tidy.core(x, row_covariates, col_covariates)
}

#' @noRd
#' @importFrom reshape2 melt
.tidy.core <- function(mats, genes, samples, genes_columns = NULL,
                       samples_columns = NULL, ...) {
  if (is.matrix(mats)) mats <- list(value = mats)
  stopifnot(is.list(mats))
  stopifnot(all(sapply(mats, is.matrix)))
  assert_named(mats, type = "unique")

  if (is.null(genes_columns)) genes_columns <- colnames(genes)
  if (is.null(samples_columns)) samples_columns <- colnames(samples)

  bad.genes <- setdiff(genes_columns, colnames(genes))
  if (length(bad.genes)) {
    stop("These columns not found in feature metadata columns: ",
         paste(bad.genes, collapse = ", "))
  }
  bad.samples <- setdiff(samples_columns, colnames(samples))
  if (length(bad.samples)) {
    stop("These columns not found in sample metadata columns: ",
         paste(bad.samples, collapse = ", "))
  }

  rnames <- rownames(mats[[1]])
  snames <- colnames(mats[[1]])
  genes$.gene_id <- rnames
  gid.col <- sapply(genes, function(xx) all(xx == rnames))
  gid.col <- colnames(genes)[which(gid.col)[1L]]
  if (gid.col != ".gene_id") genes$.gene_id <- NULL

  samples$.sample_id <- snames
  sid.col <- sapply(samples, function(xx) all(xx == snames))
  sid.col <- colnames(samples)[which(sid.col)[1L]]
  if (sid.col != ".sample_id") samples$.sample_id <- NULL

  adat.all <- lapply(names(mats), function(mname) {
    m <- mats[[mname]]
    stopifnot(all.equal(rownames(m), rnames))
    m <- reshape2::melt(m)
    m <- transform(m, Var1 = as.character(Var1), Var2 = as.character(Var2))
    colnames(m) <- c(gid.col, sid.col, mname)
    m
  })
  adat <- do.call(cbind, adat.all)
  # if there were multiple matrices, there will be multiple sample_id columns
  # so we remove those

  gcols <- unique(c(gid.col, genes_columns))
  scols <- unique(c(sid.col, samples_columns))
  adat <- adat[, !duplicated(colnames(adat))]
  out <- inner_join(adat, select(genes, !!gcols), by = gid.col)
  out <- inner_join(out, select(samples, !!scols), by = sid.col)
  as.tbl(out)
}
