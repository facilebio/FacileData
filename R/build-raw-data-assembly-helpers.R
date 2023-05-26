# ==============================================================================
# These functions are not exported
# ------------------------------------------------------------------------------
# 
# Utility functions used in wrangling raw assay data into shape for
# FacileDataSet creation.
# 
# End users will not use these functions, but they are here to whip raw
# data into shape to exercise FacileDataSet assembly and testing
# 
# Used by the `raw-assay-data-assembly.Rmd` vignette.

#' Splits a SingleCellExperiment in a list of DGELists, split by the levels
#' one of the colData columns
#' 
#' This function is called from the dataset-wrangling.Rmd vignette, and
#' relies on the SingleCellExperiment package (which is in Suggests, not
#' Imports)
#' 
#' @param x A SingleCellExperiment or SummarizedExperiment
#' @param split.by a column name in `colData(x)` that has categorical covariate
#'   (factor or character)
build_pb_split <- function(x, split.by = "cond") {
  reqpkg("SingleCellExperiment")
  split.levels <- x[[split.by]]
  stopifnot(is.character(split.levels) || is.factor(split.levels))
  
  # The current FacileDataSet constraints is that the list of datasets that
  # is used to "hydrate" each assay must all have **the same** feature space.
  y.all <- edgeR::DGEList(
    counts = SingleCellExperiment::counts(x),
    genes = as.data.frame(SummarizedExperiment::rowData(x)),
    # group was already defined in sce.all
    # group = paste0(pbx$cond, "__", pbx$cell_abbrev)
    samples = as.data.frame(SingleCellExperiment::colData(x)))
  
  # Force counts to be `integer`
  storage.mode(y.all$counts) <- "integer"
  
  y.all <- edgeR::calcNormFactors(y.all)
  
  des <- model.matrix(~ group, y.all$samples)
  keep <- edgeR::filterByExpr(y.all, design = des, min.count = 5, 
                              min.total.count = 10)
  y <- edgeR::calcNormFactors(y.all[keep,,keep.lib.sizes = FALSE])
  
  sapply(unique(split.levels), function(dsname) {
    edgeR::calcNormFactors(y[, split.levels == dsname])
  }, simplify = FALSE)
}

#' Description of available assays in this packge
#' @return a tibble of assay information
build_available_assays <- function() {
  dplyr::tribble(
    ~assay_name, ~assay_type, ~feature_type, ~storage_mode, ~description,
    "scRNAseq",  "pseudobulk", "ensgid",     "integer",     "pseudobulked scRNAseq data",
    "snRNAseq",  "pseudobulk", "ensgid",     "integer",     "pseudobulked snRNAseq data")
}

#' Returns the filepath for the raw assay data for a given assay(name)
#' @export
#' @param assay_name One of the (case insensitive) options in
#'   `available_assays()$assay_name`
#' @return file path to the serialized assaylist
build_assay_list_file_path <- function(assay_name) {
  assert_string(assay_name)
  ainfo <- build_available_assays() |> 
    mutate(name = tolower(assay_name)) |> 
    filter(.data$name == tolower(.env$assay_name))
  if (nrow(ainfo) == 0L) {
    stop("no information found for assay_name: ", assay_name)
  }
  outdir <- system.file("extdata", "assay-data", package = "FacileData")
  file.path(outdir, paste0(ainfo$assay_name, "-assay-list.qs"))
}

#' Loads an assay list object
#' @export
#' @inheritParams build_assay_list_file_path
build_assay_lists_load <- function(assay_name) {
  fn <- build_assay_list_file_path(assay_name)
  qs::qread(fn)
}
