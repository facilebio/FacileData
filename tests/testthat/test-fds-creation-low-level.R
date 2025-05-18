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

# Initialization ---------------------------------------------------------------
test_that("fds directory initialization works", {
  on.exit(cleanup_fds(fds.dir, "tfds"))  
  fds.dir <- tempfile("tes_fdsdir__")
  tfds <- fds_initialize(fds.dir, "human")
  expect_directory_exists(fds.dir)
  expect_class(tfds, "FacileDataStore")
  expect_class(tfds$con, "SQLiteConnection")
  expect_file_exists(hdf5fn(tfds), "r")
})

# Feature Space Addition -------------------------------------------------------
test_that("multiple fds feature space addition works", {
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

test_that("same feature space can't be added twice", {
  on.exit(cleanup_fds(fds.dir, "tfds"))  
  fds.dir <- tempfile("tes_fdsdir__")
  tfds <- fds_initialize(fds.dir, "human")
  
  fds_add_feature_space(
    tfds, 
    feature_info = ensembl_genes,
    feature_type = "ensgid",
    description = "example ensembl gene universe"
  )
  
  expect_error({
    fds_add_feature_space(
      tfds, 
      feature_info = ensembl_genes,
      feature_type = "ensgid",
      description = "example ensembl gene universe"
    )
  }, "already registered")
})

# Assay insertion ==============================================================
test_that("assay registration works", {
  new.assay <- tibble(
    assay = "counts",
    assay_type = "rnaseq",
    feature_type = "ensgid"
  )
  on.exit(cleanup_fds(fds.dir, "tfds"))
  
  fds.dir <- tempfile("tes_fdsdir__")
  tfds <- fds_initialize(fds.dir, "human")
  
  # obviously, no assays should registered
  ainfo.all <- assay_info(tfds)
  expect_data_frame(ainfo.all, nrows = 0L)
  
  fds_add_feature_space(
    tfds, 
    feature_info = ensembl_genes,
    feature_type = "ensgid",
    description = "example ensembl gene universe"
  )

  fds_register_assay(
    tfds,
    assay_name = new.assay$assay,
    assay_type = new.assay$assay_type,
    feature_type = new.assay$feature_type,
    storage_mode = "integer"
  )
  
  # check that the assay group was created in the hdf5 file
  fid <- rhdf5::H5Fopen(hdf5fn(tfds))
  info <- tryCatch({
    rhdf5::H5Gget_info_by_name(fid, "/assay/counts")
  }, error = function(x) NULL)
  expect_list(info, names = "unique")
  
  # check db tables setup
  ainfo <- assay_info(tfds)
  expect_data_frame(ainfo, nrows = 1L)
  expect_equal(
    select(ainfo, all_of(colnames(new.assay))),
    new.assay, 
    check.attributes = FALSE)
  
  afeatures <- features(tfds, assay = new.assay$assay)
  expect_data_frame(afeatures, nrows = nrow(ensembl_genes))
  
  # ensure entire 1:nrow() sequential space covered by hdf5_index column
  expect_set_equal(afeatures$hdf5_index, 1:nrow(afeatures))
  expect_equal(
    select(afeatures, feature_id, name, meta),
    select(ensembl_genes, feature_id, name, meta),
    check.attributes = FALSE
  )
})

test_that("duplicate assay registration throws an error", {
  new.assay <- tibble(
    assay = "counts",
    assay_type = "rnaseq",
    feature_type = "ensgid"
  )
  on.exit(cleanup_fds(fds.dir, "tfds"))
  
  fds.dir <- tempfile("tes_fdsdir__")
  tfds <- fds_initialize(fds.dir, "human")
  
  fds_add_feature_space(
    tfds, 
    feature_info = ensembl_genes,
    feature_type = "ensgid",
    description = "example ensembl gene universe"
  )
  
  fds_register_assay(
    tfds,
    assay_name = new.assay$assay,
    assay_type = new.assay$assay_type,
    feature_type = new.assay$feature_type,
    storage_mode = "integer"
  )
  
  expect_error({
    fds_register_assay(
      tfds,
      assay_name = new.assay$assay,
      assay_type = new.assay$assay_type,
      feature_type = new.assay$feature_type,
      storage_mode = "integer"
    )
  }, "already exists")
})

test_that("fds_add_assay_data maintains fidelity of data", {
  new.assay <- tibble(
    assay = "counts",
    assay_type = "rnaseq",
    feature_type = "ensgid"
  )
  se <- se.list[[1]]
  se.counts <- SummarizedExperiment::assay(se, "counts")
  
  on.exit(cleanup_fds(fds.dir, "tfds"))
  fds.dir <- tempfile("tes_fdsdir__")
  tfds <- fds_initialize(fds.dir, "human")
  fspace <- fds_add_feature_space(
    tfds, 
    feature_info = ensembl_genes,
    feature_type = "ensgid",
    description = "example ensembl gene universe"
  )
  rassay <- fds_register_assay(
    tfds,
    assay_name = new.assay$assay,
    assay_type = new.assay$assay_type,
    feature_type = new.assay$feature_type,
    storage_mode = "integer"
  )
  
  ad <- fds_add_assay_data(
    tfds,
    se.counts,
    assay_name = new.assay$assay,
    dataset_name = se$dataset[1]
  )
  
  amatrix <- fetch_assay_data(tfds, as.matrix = TRUE)
  colnames(amatrix) <- sub(".*?__", "", colnames(amatrix))
  
  expect_equal(dim(amatrix), dim(se))  
  expect_setequal(rownames(amatrix), rownames(se))
  expect_setequal(colnames(amatrix), colnames(se))
  expect_equal(
    amatrix[rownames(se), colnames(se)],
    se.counts,
    check.attributes = FALSE)
})
