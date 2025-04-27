context("Feature Descriptors")

if (!exists("afds")) afds <- an_fds()
.fids.good <- c("ENSG00000000460", "ENSG00000001167", "ENSG00000000005")
.fids.bad <- c("ESNG002", "ENS00004")

test_that("undefined assay_name uses default_assay", {
  # If we don't specify an assay_name, the default_assay will be returned
  fd <- create_assay_feature_descriptor(afds, .fids.good)
  expect_setequal(fd$feature_id, .fids.good)
  expect_equal(unique(fd$assay), default_assay(afds))
})

test_that("duplicated feature_ids/assay_name are removed", {
  # Duplicated 
  fd <- create_assay_feature_descriptor(afds, rep(.fids.good, 2))
  expect_equal(nrow(fd), length(.fids.good))
  
    
  fids.aname <- tibble(
    feature_id = rep(.fids.good, 2),
    assay = rep(assay_names(afds)[1:2], each = length(.fids.good)))

  fd <- create_assay_feature_descriptor(afds, fids.aname)
  expect_equal(nrow(fd), length(.fids.good) * 2)         
  expect_setequal(
    paste(fd$feature_id, fd$assay),
    paste(fids.aname$feature_id, fids.aname$assay))
})

test_that("alternate assay_name supporrted", {
  alt.assay <- setdiff(assay_names(afds), default_assay(afds))[1]
  fd <- create_assay_feature_descriptor(
    afds,
    .fids.good,
    assay_name = alt.assay
  )
  expect_setequal(fd$feature_id, .fids.good)
  expect_equal(unique(fd$assay), alt.assay)
  expect_false(unique(fd$assay) == default_assay(afds))
})

test_that("no defined feature_id returns all features for assay", {
  fids.default <- create_assay_feature_descriptor(afds)
  expected.default <- features(afds)
  expect_setequal(fids.default$feature_id, expected.default$feature_id)
  expect_equal(unique(fids.default$assay), default_assay(afds))

  alt.assay <- setdiff(assay_names(afds), default_assay(afds))[1]
  expect_equal(alt.assay, assay_names(afds)[2])  
  fids.alt <- create_assay_feature_descriptor(afds, assay_name = alt.assay)
  expected.alt <- features(afds, assay_name = alt.assay)
  expect_setequal(fids.alt$feature_id, expected.alt$feature_id)
})
