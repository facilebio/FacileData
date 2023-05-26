# A set of tests around building multi-modal faciledataset objects.
# these are meant to more thoroughly test different aspects of dataset assembly
# and use the "well curated" set of 'assay data list objects' assembled in the
# `raw-assay-data-assembly.Rmd` vignette.

test_that("A list of assay data objects assemble into single-assay FacileDataSet", {
  name <- "TestSingleAssayFacileDataSet"
  ainfo <- build_available_assays() |> 
    filter(assay_name == "scRNAseq")
  adat <- build_assay_lists_load(ainfo$assay_name)
  
  path <- tempfile(paste0(name, "_"))
  fds <- as.FacileDataSet(
    adat,
    path = path,
    dataset_name = name,
    assay_name = ainfo$assay_name,
    assay_type = ainfo$assay_type,
    assay_description = ainfo$description,
    organism = "Homo sapiens")
  
  # check feature space --------------------------------------------------------
  
  # Currently we assume that all of the features of the individuals datasets
  # that were passesd into the `adat` list to create the dataset had the same
  # feature space, so let's assert that's still true here.
  for (features.y in lapply(adat, "[[", "genes")) {
    expect_equal(rownames(features.y), rownames(adat[[1]]), info = "just checking")
  }
    
  features.fds <- features(fds)
  expect_equal(nrow(features.fds), nrow(features.y))
  expect_set_equal(features.fds$feature_id, features.y$feature_id)
  # We'll take advantage of the fact that rownames(features.y) can be used
  # as an index into the rows
  fy <- features.y[features.fds$feature_id,]
  check.cols <- c("feature_id", "feature_type", "name", "meta")
  expect_equal(features.fds[, check.cols], fy[, check.cols], 
               check.attributes = FALSE)
  
  # check assay data was transferred faithfully --------------------------------
  for (dsname in names(adat)) {
    y <- adat[[dsname]]
    # let's grab some random samples
    sidx <- sample(ncol(y), min(5, ncol(y)))
    sids <- y$samples$sample_id[sidx]
    # some random features
    ridx <- sample(nrow(y), min(5, nrow(y)))
    
    expected <- y$counts[ridx, sidx]
    
    # facilize
    dsamples <- tibble(dataset = dsname, sample_id = sids)
    dfeatures <- y$genes[ridx,]
    # sanity check
    expect_equal(colnames(expected), dsamples$sample_id)
    expect_equal(rownames(expected), dfeatures$feature_id)
    
    # now check the assay retrieval
    fdat <- fetch_assay_data(fds, dfeatures, dsamples, as.matrix = TRUE)
    # let's strim the dataset__ column name prefix
    colnames(fdat) <- sub(paste0(dsname, "__"), "", colnames(fdat))
    
    expect_set_equal(rownames(fdat), rownames(expected))
    expect_set_equal(colnames(fdat), colnames(expected))
    expect_equal(fdat[rownames(expected), colnames(expected)], expected)
  }
})
    