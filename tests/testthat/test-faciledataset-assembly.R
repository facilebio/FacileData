# A set of tests around building multi-modal faciledataset objects.
# these are meant to more thoroughly test different aspects of dataset assembly
# and use the "well curated" set of 'assay data list objects' assembled in the
# `raw-assay-data-assembly.Rmd` vignette.

test_that("List of assay data objects assemble into single-assay FacileDataSet", {
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

test_that("Supports addition of second assay to original feature new features", {
  name <- "TestMultiAssayFacileDataSet"
  path <- tempfile(paste0(name, "_"))

  assays <- build_available_assays()
  ai1 <- assays[1,]
  ad1 <- build_assay_lists_load(ai1$assay_name)
  fids1 <- lapply(ad1, \(x) rownames(x)) |> unlist() |> unique()
  
  fds <- as.FacileDataSet(
    ad1,
    path = path,
    dataset_name = name,
    assay_name = ai1$assay_name,
    assay_type = ai1$assay_type,
    assay_description = ai1$description,
    organism = "Homo sapiens")
  f.original <- features(fds, feature_type = ai1$feature_type)
  
  # We've already tested this, but let's make sure the feature space in the
  # FacileDataSet matches the feature space of the data we sent in already
  expect_set_equal(f.original$feature_id, fids1)
  
  # Now add the second assay ...
  ad2 <- build_assay_lists_load(ai2$assay_name)
  # ... which is also over the `"ensgid"` feature space.
  ai2 <- assays[2,]
  expect_equal(ai1$feature_type, ai2$feature_type)

  # ... but the overlap of feature_id's from first dataset to this one
  # is not perfect.
  fids2 <- lapply(ad2, \(x) rownames(x)) |> unlist() |> unique()
  
  # 1. There is a large number of features that are shaerd:
  expect_gt(length(intersect(fids1, fids2)), 15000)
  
  # 2. There are some features in the new assay that aren't in the old one
  f2.new <- ad2[[1]]$genes |> 
    anti_join(f.original, by = "feature_id") |> 
    as_tibble()
  expect_gt(nrow(f2.new), 5000)

  # Now we add the assay data, and the new features should be added to the
  # database
  addFacileAssaySet(
    fds,
    ad2,
    facile_assay_name = ai2$assay_name,
    facile_assay_type = ai2$assay_type,
    facile_feature_type = ai2$feature_type,
    facile_assay_description = ai2$description,
    facile_feature_info = ad2[[1]]$genes,
    storage_mode = ai2$storage_mode,
    assay_name = ai2$assay_name)
  
  f.universe <- features(fds, feature_type = "ensgid")
  
  expect_gt(nrow(f.universe), nrow(f.original))
  expect_subset(f.original$feature_id, f.universe$feature_id)
  expect_subset(f2.new$feature_id, f.universe$feature_id)
  expect_set_equal(c(f.original$feature_id, f2.new$feature_id), f.universe$feature_id)
})
