#' Assembles an example facile dataset to play with
#'
#' This combines the airway and parathyroidSE RNA-seq datasets into a single
#' FacileDataSet.
#'
#' The code here is extracted from the `FacileDataSet-assembly` vignette. Please
#' read that for some of the why's and how's of the decisions made here when
#' assembling datasets.
#'
#' @export
#' @param directory The name of the parent directory to hold the dataset
#' @param name A subdirectory within `directory` will be created using this
#'   name.
#' @return The FacileDataSet object itself.
#' @examples
#' \dontest{
#' afds <- assemble_example_dataset()
#' }
assemble_example_dataset <- function(directory = tempdir(),
                                     name = "ExampleRnaFacileDataSet") {
  assert_directory_exists(directory, "w")
  full.path <- file.path(directory, name)
  if (file.exists(full.path)) {
    stop("The output directory already exists, remove it to recreate the ",
         "dataset:\n  ", full.path)
  }
  message("Assembling dataset into: ", full.path)

  reqpkg("SummarizedExperiment")

  # Load Data ..................................................................
  dat.env <- new.env()
  utils::data("airway", package = "airway", envir = dat.env)
  utils::data("parathyroidGenesSE", package = "parathyroidSE", envir = dat.env)

  # Munge colData ..............................................................
  se.airway <- dat.env[["airway"]]
  cd.airway <- SummarizedExperiment::colData(se.airway) |>
    as.data.frame() |>
    transmute(
      sample_type = "cell_line",
      cell_line = cell,
      treatment = ifelse(dex == "untrt", "control", "dex"))
  rownames(cd.airway) <- colnames(se.airway)

  se.parathyroid <- dat.env[["parathyroidGenesSE"]]
  cd.parathyroid <- SummarizedExperiment::colData(se.parathyroid) |>
    as.data.frame() |>
    transmute(
      sample_type = "primary",
      subject_id = paste0("patient_", patient),
      treatment = tolower(as.character(treatment)),
      time = paste0("hrs", sub("h$", "", time)))
  rownames(cd.parathyroid) <- colnames(se.parathyroid)

  # Munge rowData ..............................................................
  mart.info <- local({
    fn <- system.file("extdata", "ensembl-v75-gene-info.csv.gz",
                      package = "FacileData")
    con <- gzfile(fn, "rt")
    on.exit(close.connection(con))
    read.csv(con, stringsAsFactors = FALSE)
  })

  shared.ids <- intersect(rownames(se.airway), rownames(se.parathyroid))
  gene.info <- mart.info |>
    transmute(feature_id = ensembl_gene_id,
              feature_type = "ensgid",
              name = hgnc_symbol,
              meta = gene_biotype,
              source = "Ensembl_v75") |>
    filter(feature_id %in% shared.ids) |>
    distinct(feature_id, .keep_all = TRUE)
  rownames(gene.info) <- gene.info[["feature_id"]]

  y.airway <- edgeR::DGEList(
    counts = SummarizedExperiment::assay(se.airway)[rownames(gene.info), ],
    samples = cd.airway,
    genes = gene.info)

  y.para <- edgeR::DGEList(
    counts = SummarizedExperiment::assay(se.parathyroid)[rownames(gene.info), ],
    samples = cd.parathyroid,
    genes = gene.info)

  dat.all <- list(airway = y.airway, parathyroid = y.para)

  xfds <- as.FacileDataSet(
    dat.all,
    path = full.path,
    dataset_name = name,
    assay_name = "gene_counts",
    assay_description = "Gene counts provided by Bioconductor data packages",
    assay_type = "rnaseq",
    organism = "Homo sapiens")

  xfds
}
