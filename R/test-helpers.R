# multimodal faciledataset functions ===========================================
# 
# The multimodal faciledataset was assembled by pseudo-bulking the scRNAseq
# and snRNAseq data from teh KPMP into their own assays.
# 
# The dataset,sample_id pairs look like:
#    dataset sample_id     
#    <chr>   <chr>         
#  1 HCKD    PT.31_10000   
#  2 HCKD    DCT.31_10000  
#  3 DKD     PT.31_10001   
#  4 DKD     DCT.29_10012  
#
# dataset = "condition" (healthy or otherwise)
# sample_id = celltype<>donor_id

#' Retrieves a multi-modal FacileDataSet based on KPMP data.
#' 
#' We use this dataset for examples and unit testing purposes.
#' 
#' @export
#' @return a FacileDataSet
an_fds <- function() {
  fn <- system.file("extdata", "testFacileDataSet", package = "FacileData")
  FacileDataSet(fn)
}

#' A set of features that show variable expression across some_celltypes()
#' 
#' These features had strong ANOVA hits across [some_celltypes()]
#' 
#' @export
some_features <- function() {
  out <- tibble::tribble(
    ~feature_id,        ~name,
    "ENSG00000119888",  "EPCAM",
    "ENSG00000149564",  "ESAM",
    "ENSG00000147676",  "MAL2",
    "ENSG00000128567",  "PODXL")
  mutate(out, assay = "scrnaseq")
}

#' A set of cell types that exhibit strong expression variability.
#' 
#' The [testing_features()] show cell type expression patterns across these
#' celltypes
#' 
#' @export
some_celltypes <- function() {
  c("B", "C-PC", "CNT", "DCT", "PT", "TAL", "EC", "MPC", "POD")
}

#' Retrieve a set of samples for testing
#' 
#' @export
#' @param sparse If `FALSE` (default), a set of samples are returned that have
#'   assay data from both assays in the dataset, otherwise some samples will
#'   be missing scRNAseq data.
some_samples <- function(sparse = FALSE) {
  fds. <- an_fds()
  if (sparse) {
    # Manual inspection of the samples an an_fds() shows that the selection
    # criteria below provides a list of samples that all have data from the
    # scrnaseq assay
    out <- filter_samples(fds., cell_abbrev %in% some_celltypes())
  } else {
    # If we exclude donor 29-10008, the PT and IMM celltypes have assay data
    # for all samples from the scRNAseq and snRNAseq assay
    out <- filter_samples(fds.,
                          cell_abbrev %in% c("PT", "IMM"),
                          donor_id != "29-10008")
  }
  out
}

# ------------------------------------------------------------------------------

#' Creata a FacileDataSet from the individual assay lists
#' 
#' @export
#' @param name The name of the datasets. Defaults to
#'   `"TestMultiModalFacileDataSet"`. Users may create different versions to
#'   test (ie. different assays to include, etc) and may want a more
#'   recognizable name for it.
#' @param path Local path to create the FacileDataSet directory. If `NULL`
#'   (default) a temporary directory will be created via`tempfile()` and will
#'   include `name` as the directory's prefix.
#' @param assays the nums of the assays to include. If `NULL` (default), all
#'   assays will be included.
assembleFacileDataSet <- function(name = "TestMultiModalFacileDataSet",
                                  path = NULL, assays = NULL) {
  if (FALSE) {
    name <- "TestMultiModalFacileDataSet"
    path <- NULL
    assays <- NULL
  }
  assert_string(name)
  if (is.null(path)) path <- tempfile(paste0(name, "_"))
  assert_directory(dirname(path), "w")
  
  all_assays <- build_available_assays()
  if (is.null(assays)) assays <- all_assays$assay_name
  assays <- assert_subset(tolower(assays), tolower(all_assays$assay_name))
  
  adata <- all_assays |> 
    mutate(assay_name = tolower(assay_name)) |> 
    semi_join(tibble(assay_name = assays), by = "assay_name")
  
  ainfo <- adata[1,]
  adat <- build_assay_lists_load(ainfo$assay_name)
  
  fds <- as.FacileDataSet(
    adat,
    path = path,
    dataset_name = name,
    assay_name = ainfo$assay_name,
    assay_type = ainfo$assay_type,
    assay_description = ainfo$description,
    organism = "Homo sapiens")
  
  if (nrow(adata) > 1L) {
    for (i in 2:nrow(adata)) {
      ainfo <- adata[i,]
      adat <- build_assay_lists_load(ainfo$assay_name)
      
      if (is(adat[[1L]], "DGEList")) {
        assay_name <- "counts"
      } else if (is(adat[[1L]], "EList")) {
        assay_name <- "E"
      } else {
        stop("No handler for class yet: ", class(adat[[1]])[1L])
      }
      message("Adding assay: ", ainfo$assay_name)
      addFacileAssaySet(
        fds,
        adat,
        facile_assay_name = ainfo$assay_name,
        facile_assay_type = ainfo$assay_type,
        facile_feature_type = ainfo$feature_type,
        facile_assay_description = ainfo$description,
        facile_feature_info = adat[[1]]$genes,
        storage_mode = ainfo$storage_mode,
        assay_name = assay_name)
    }
  }
  
  fds
}

# Old Helpers ==================================================================
# These were based on manipulation of a precanned FacileDataSet based on TCGA
# data that came prepackagd with this package.
# 
# There was no code or data made available to rebuild this dataset. It was
# useful for testing, but we want to phase it out.

#' Retrieves an example FacileDataSet
#'
#' A subset of the TCGA data from the BLCA and COAD indications is provided
#' as a FacileDataSet.
#'
#' @export
exampleFacileDataSet <- function() {
  fn <- system.file('extdata', 'exampleFacileDataSet', package='FacileData')
  FacileDataSet(fn)
}

#' Fetches exemplar data for unit testing
#'
#' @export
#' @rdname test-helpers
example_sample_covariates <- function() {
  pdat <- system.file("testdata", "test-sample-covariates.rds",
                      package = "FacileData")
  readRDS(pdat)
}

#' @export
#' @rdname test-helpers
#' @param file.path If `TRUE`, returns the path to the yaml file, otherwise
#'   returns the list-of-list meta definition.
#' @return Either the list-of-list meta definition, or path to the `meta.yaml`
#'   file where these are defined.
example_meta <- function(file.path=FALSE) {
  out <- system.file("testdata", "expected-meta.yaml",
                     package = "FacileData")
  if (!isTRUE(file.path)) {
    out <- yaml::read_yaml(out)
  }
  out
}

#' @export
#' @importFrom yaml read_yaml
#' @rdname test-helpers
#' @return the list-of-list definitions for the example `pData` returned from
#'   [example_sample_covariates()]
example_sample_covariate_definitions <- function() {
  out <- example_meta(file.path=FALSE)
  out$sample_covariate
}
