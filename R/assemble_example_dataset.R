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
assemble_example_dataset <- function(directory = tempdir(),
                                     name = "ExampleRnaFacileDataSet") {
  assert_directory_exists(directory, "w")
  full.path <- file.path(directory, name)
  if (file.exists(full.path)) {
    stop("The output directory already exists, remove it to recreate the ",
         "dataset:\n  ", full.path)
  }
  message("Assembling dataset into: ", full.path)

  ns <- tryCatch(loadNamespace("SummarizedExperiment"), error = function(e) NULL)
  if (is.null(ns)) stop("SummarizedExperiment required")
  ns4 <- tryCatch(loadNamespace("S4Vectors"), error = function(e) NULL)
  if (is.null(ns4)) stop("S4Vectors required")

  # Load Data ..................................................................
  dat.env <- new.env()
  tryCatch({
    data("airway", package = "airway", envir = dat.env)
  }, error = function(e) stop("The airway package is required"))
  tryCatch({
    data("parathyroidGenesSE", package = "parathyroidSE", envir = dat.env)
  }, error = function(e) stop("The parathyroidSE package is required"))

  # Munge colData ..............................................................
  se.airway <- dat.env[["airway"]]
  cd.airway <- local({
    cd <- ns$colData(dat.env[["airway"]]) %>%
      as.data.frame() %>%
      transmute(
        sample_type = "cell_line",
        cell_line = cell,
        treatment = ifelse(dex == "untrt", "control", "dex")) %>%
      ns4$DataFrame()
    rownames(cd) <- colnames(se.airway)
    cd
  })
  se.airway <- ns$`colData<-`(se.airway, value = cd.airway)

  se.parathyroid <- dat.env[["parathyroidGenesSE"]]
  cd.parathyroid <- local({
    cd <- ns$colData(dat.env[["parathyroidGenesSE"]]) %>%
      as.data.frame() %>%
      transmute(
        sample_type = "primary",
        subject_id = paste0("patient_", patient),
        treatment = tolower(as.character(treatment)),
        time = paste0("hrs", sub("h$", "", time))) %>%
      ns4$DataFrame()
    rownames(cd) <- colnames(se.parathyroid)
    cd
  })
  se.parathyroid <- ns$`colData<-`(se.parathyroid, value = cd.parathyroid)

  # Munge rowData ..............................................................
  mart.info <- local({
    fn <- system.file("extdata", "ensembl-v75-gene-info.csv.gz",
                      package = "FacileData")
    read.csv(gzfile(fn, "rt"), stringsAsFactors = FALSE)
  })

  shared.ids <- intersect(rownames(se.airway), rownames(se.parathyroid))
  gene.info <- mart.info %>%
    transmute(feature_id = ensembl_gene_id,
              feature_type = "ensgid",
              name = hgnc_symbol,
              meta = gene_biotype,
              source = "Ensembl_v75") %>%
    filter(feature_id %in% shared.ids) %>%
    distinct(feature_id, .keep_all = TRUE)
    # ns4$DataFrame()
  rownames(gene.info) <- gene.info[["feature_id"]]

  # I assemble these into DGELists because I can't figure out how to get
  # SummarizedExperiment subsetting working without using the loadedNamespace
  # mojo ... I'm dying on the inside here.
  #
  # Obviously I'm doing something wrong, but ... damn, y0 ... damn.

  # se.subfn <- selectMethod("[", c("SummarizedExperiment", "ANY", "ANY"))

  # se.airway <- se.airway[rownames(gene.info),]
  # se.airway <- se.subfn(se.airway, rownames(gene.info))
  # se.airway <- ns$`rowData<-`(se.airway, value = gene.info)
  #
  # # se.parathyroid <- se.parathyroid[rownames(gene.info),]
  # se.parathyroid <- se.subfn(se.parathyroid, rownames(gene.info))
  # se.parathyroid <- ns$`rowData<-`(se.parathyroid, value = gene.info)
  #
  # dat.all <- list(airway = se.airway, parathyroid = se.parathyroid)

  # gene.info <- ns4$as.data.frame(gene.info)

  y.airway <- edgeR::DGEList(
    counts = ns$assay(se.airway)[rownames(gene.info), ],
    samples = ns4$as.data.frame.DataTable(ns$colData(se.airway)),
    genes = gene.info)

  y.para <- edgeR::DGEList(
    counts = ns$assay(se.parathyroid)[rownames(gene.info), ],
    samples = ns4$as.data.frame.DataTable(ns$colData(se.parathyroid)),
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
