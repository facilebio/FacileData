#' Create a Facile___DataSet Package from a list of SummarizedExperiments
#'
#' Accepts a list of datasets (SummarizedExperiments) and metadata and
#' builds/installs a Facile___DataSet based on this.
#'
#' @md
#' @param data_list named `list` of data `SummarizedExperiment` objects to be converted into a `Facile___DataSet` package.
#' @param slug name to use in the name `Facile___DataSet` package.
#' @param version version to label the new `Facile___DataSet` package.
#' @param parent_path directory in which the package will be created. Must already exist.
#' @param covariates covariates to extract from the listed data structures. Must exist in each dataset.
#' @param cov_metadata nested `list` of covariate metadata. Must contain at least label, description, type and class for each covariate.
#' @param data_metadata named `list`` of metadata for the datasets. Must contain at least url and description for each dataset.
#' @param assay_name label to be used within FacileExplorer. Default "rnaseq".
#' @param organism label to be used within FacileExplorer. Default "Homo Sapiens".
#'
#' @importFrom S4Vectors mcols metadata
#' @importFrom glue glue
#' @importFrom SummarizedExperiment colData
#' @importFrom devtools create document install build
#' @importFrom desc desc_set_version
#'
#' @return Builds and installs a `Facile___DataSet` package, leaving a source tarball in `parent_path` directory.
#' @export
#'
#' @examples \dontrun{
#' ngs114 <- ep.project.from.id("PRJ0013166", standard.checks = FALSE)
#' ngs171 <- ep.project.from.id("PRJ0013155", standard.checks = FALSE)
#'
#' data_list <- list(
#'   NGS114 = ep.SummarizedExperiment(ngs114, attach.annot = TRUE),
#'   NGS171 = ep.SummarizedExperiment(ngs171, attach.annot = TRUE)
#' )
#'
#' cov_metadata <- list(
#'   primary_tissue = list(label = "Primary Tissue",
#'                         description = "Primary tissue source",
#'                         type = "tumor_classification",
#'                         class = "categorical"
#'   ),
#'   gender = list(label = "Gender/Sex",
#'                 description = "In the 'ratio between chrX:chrY' sense.",
#'                 type = "clinical",
#'                 class = "categorical"
#'   ),
#'   diagnosis = list(label = "Tissue Diagnosis",
#'                    description = "Tissue metaclass oncology",
#'                    type = "tumor_classification",
#'                    class = "categorical"
#'   ),
#'   ethnicity = list(label = "Ethnicity",
#'                    description = "Patient ethnicity",
#'                    type = "clinical",
#'                    class = "categorical"
#'   ),
#'   tissue = list(label = "Tissue Group",
#'                 description = "Rollup of tissue type to defined vocab",
#'                 type = "tumor_classification",
#'                 class = "categorical"
#'   )
#' )
#'
#' data_metadata <- list(
#'   NGS114 = list(url = "http://gene.com", description = "This is NGS114"),
#'   NGS171 = list(url = "http://gene.com", description = "This is NGS171")
#' )
#'
#' create_FDS_pkg(data_list = data_list,
#'                slug = "GCell",
#'                version = "0.0.1",
#'                parent_path = "~/FacileVerse",
#'                covariates = c("PRIMARY_TISSUE", "GENDER",
#'                               "TISSUE_DIAGNOSIS", "ETHNICITY",
#'                               "TISSUE_METACLASS_ONCOLOGY"),
#'                cov_metadata = cov_metadata,
#'                data_metadata = data_metadata)
#'
#' }
create_FDS_pkg <- function(data_list = NULL,
                           slug = NULL,
                           version = "0.0.1",
                           parent_path = ".",
                           covariates = NULL,
                           cov_metadata = NULL,
                           data_metadata = NULL,
                           assay_name = "rnaseq",
                           source_assay = "counts",
                           organism = "Homo sapiens") {

  if (is.null(slug)) stop("Please provide a slug to use in the new dataset name: Facile<slug>DataSet")
  FDS_name <- glue::glue("Facile{slug}DataSet")
  FDS_version <- version

  ## assertions
  if (is.null(covariates)) stop("Please specify covariates to extract")
  if (is.null(cov_metadata)) stop("Please provide metadata for each selected covariate")
  if (!identical(length(covariates), length(cov_metadata))) stop("cov_metadata must have the same length as covariates")
  if (is.null(data_metadata)) stop("Please provide metadata for each dataset")
  if (!identical(length(data_list), length(data_metadata))) stop("data_metadata must have an element for each dataset")
  if (!identical(names(data_list), names(data_metadata))) stop("data_list and data_metadata must have identical names")
  if (!file.exists(parent_path)) stop("Directory to contain new FDS directory must already exist")

  ## examples of type: data_batch, clinical, tumor_classification, response, IHC, mutation, safety, conmeds
  ## examples of class: numerical, categorical, right_censored
  ## (**can't create new type or class**)
  curated_data_list <- lapply(seq_along(data_list), function(i, d) {

    ## work with the current list element only
    x <- d[[i]]

    ## assertions
    if (!inherits(x, "SummarizedExperiment")) stop("Requires SummarizedExperiment objects in data_list")
    if (length(colnames(mcols(x))) != 5L) stop("Requires exactly 5 columns in rowData")
    if (any(!(covariates %in% names(cov_metadata)))) stop("All covariates must be present in colData for all datasets")

    ## prepend feature_id with "GeneID:"
    # if (!grepl("GeneID", rownames(x)[1])) {
    #   rownames(x) = paste0("GeneID:", rownames(x))
    # }

    new_colnames <- c("aliases", "effective_length", "feature_type", "name", "meta")
    if (!identical(colnames(mcols(x)), new_colnames)) {
      if (i == 1) { # quiet on subsequent elements
        message("* Renaming rowData columns:")
        invisible(sapply(paste(colnames(mcols(x)), new_colnames, sep = " --> "), message))
      }
      colnames(mcols(x)) = new_colnames
    }

    ## source is a new column
    ## feature_type is overwritten
    mcols(x)$source = "IGIS"
    mcols(x)$feature_type = "entrez"
    mcols(x)$feature_id = rownames(x)

    ## ensure that all cov_metadata entries contain minimal metadata entries
    required_metadata <- c("label", "description", "type", "class")
    all(
      sapply(
        lapply(cov_metadata, names),
        function(x) {
          all(required_metadata %in% x)
        }
      )
    )

    colData(x) <- colData(x)[, covariates]
    colnames(colData(x)) <- names(cov_metadata)
    metadata(colData(x)) <- cov_metadata

    ## ensure that all data_metadata entries contain minimal metadata entries
    required_data_metadata <- c("url", "description")
    all(
      sapply(
        lapply(data_metadata, names),
        function(x) {
          all(required_data_metadata %in% x)
        }
      )
    )
    metadata(x) <- data_metadata[[names(d)[i]]]

    return(x)

  }, d = data_list)

  ## reset names
  names(curated_data_list) <- names(data_list)

  ## create new directory to store package
  DIR <- file.path(parent_path, FDS_name)
  unlink(DIR, recursive = TRUE, force = TRUE)
  devtools::create(DIR) ## create a package
  ## start an installation directory
  ## the parent directory of the
  ## FacileDataSet must already exist.
  dir.create(file.path(DIR, "inst", "extdata"), recursive = TRUE)

  ## create FDS
  message(paste0("* Creating ", FDS_name, " object"))
  message("* This can take some time.")
  as.FacileDataSet(curated_data_list,
                   path = file.path(DIR, "inst", "extdata", FDS_name),
                   assay_name = assay_name,
                   assay_type = "rnaseq", # required
                   source_assay = source_assay,
                   dataset_name = slug,
                   organism = organism
  )

  ## create accessor function dynamically
  writeLines(
    glue::glue("##' A connection to the {FDS_name}",
               "#'",
               "#' @export",
               "#' @return A \\code{{FDS_name}} object",
               "{FDS_name} <- function(path, cache_size=80000,",
               "                               db.loc = c('reference', 'temporary', 'memory')) {{",
               "  if (missing(path) || is.null(path)) {{",
               "    path <- system.file('extdata', '{FDS_name}', package='{FDS_name}')",
               "  }}",
               "  db.loc <- match.arg(db.loc)",
               "  FacileDataSet(path, cache_size = cache_size, db.loc = db.loc)",
               "}}", .sep = "\n"),
    file.path(DIR, "R", glue::glue("{FDS_name}.R"))
  )

  ## custom-annotation directory would be removed if empty
  writeLines("", file.path(DIR, "inst", "extdata", FDS_name, "custom-annotation", "README.md"))
  ## update package version
  desc::desc_set_version(FDS_version, file = file.path(DIR, "DESCRIPTION"))
  ## document, install, and build source tar.gz
  devtools::document(DIR)
  devtools::install(DIR)
  devtools::build(pkg = DIR, path = parent_path)

}
