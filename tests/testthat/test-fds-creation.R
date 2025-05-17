'
se.list <- "list of SummarizedExperiments to create a FDS from"

xfds <- fds_initialize(path, "human")
fds_register_feature_space(xfds, ensembl_genes, feature_type = "ensgid")

fds_register_assay(
  xfds,
  assay_name   = "counts",
  assay_type   = "rnaseq",
  feature_type = "ensgid",
  storage_mode = "integer")

fds_register_assay(
  xfds,
  assay_name   = "abundance",
  assay_type   = "tpm",
  feature_type = "ensgid",
  storage_mode = "integer")

for (se in se.list) {
  fds_add_assay_data(
    xfds,
    assay(se, "counts"),
    dataset_name = se$dataset[1],
    assay_name = "counts",
  )
  fds_add_assay_data(
    xfds,
    assay(se, "abundance"),
    dataset_name = se$dataset[1],
    assay_name = "abundance",
  )
  fds_add_sample_data(
    xfds,
    colData(se),
    dataset_name = se$dataset[1]
  )
}
'

if (!exists("se.list")) {
  se.list <- an_se_list()
  ensembl_genes <- SummarizedExperiment::rowData(se.list[[1L]]) |> 
    as_tibble() |> 
    transmute(feature_id, name, meta, feature_type, source = "example fds") |> 
    arrange(feature_id)
}

cleanup_fds <- function(fds.dir, fds.name = "tfds", ..., verbose = FALSE) {
  if (exists(fds.name, parent.frame())) {
    if (verbose) message("Cleaning up temp fds: ", fds.name)
    xx <- get(fds.name, parent.frame())
    if (is(xx, "FacileDataStore")) {
      if (is(xx$con, "SQLiteConnection")) RSQLite::dbDisconnect(xx$con)
    }
  }
  if (dir.exists(fds.dir)) unlink(fds.dir, recursive = TRUE)
}

test_that("fds directory initialization works", {
  on.exit(cleanup_fds(fds.dir, "tfds"))  
  fds.dir <- tempfile("tes_fdsdir__")
  tfds <- fds_initialize(fds.dir, "human")
  expect_directory_exists(fds.dir)
  expect_class(tfds, "FacileDataStore")
  expect_class(tfds$con, "SQLiteConnection")
})

test_that("fds feature space addition works", {
  on.exit(cleanup_fds(fds.dir, "tfds"))  
  fds.dir <- tempfile("tes_fdsdir__")
  tfds <- fds_initialize(fds.dir, "human")
  
  fds_add_feature_space(
    tfds, 
    feature_info = ensembl_genes,
    feature_type = "ensgid",
    description = "example ensembl gene universe"
  )
  
  features.fds <- features(tfds, feature_type = "ensgid") |> 
    arrange(feature_id) |> 
    select(all_of(colnames(ensembl_genes)))
  expect_equal(features.fds, ensembl_genes, check.attributes = FALSE)
  
  # override feature_type
  fds_add_feature_space(
    tfds, 
    feature_info = ensembl_genes,
    feature_type = "ensgid2",
    description = "example ensembl gene universe2"
  )
  
  expect_setequal(feature_types(tfds), c("ensgid", "ensgid2"))
  features2.fds <- features(tfds, feature_type = "ensgid2") |> 
    arrange(feature_id) |> 
    select(all_of(colnames(ensembl_genes)))
  expect_equal(unique(features2.fds$feature_type), "ensgid2")
  expect_equal(
    select(features2.fds, -feature_type),
    select(ensembl_genes, -feature_type),
    ensembl_genes, check.attributes = FALSE)
})
