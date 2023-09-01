if (!exists("ssamples")) ssamples <- FacileData::some_samples(sparse = TRUE)
if (!exists("fsamples")) fsamples <- FacileData::some_samples(sparse = FALSE)

test_that("has_assay decoreates samples with assay status", {
  prefix <- "has_"
  anames <- assay_names(ssamples)
  expected_colnames <- paste0(prefix, anames)
  
  has <- has_assay(ssamples, prefix = prefix)
  expect_equal(nrow(has), nrow(ssamples))
  expect_subset(expected_colnames, colnames(has))
  
  for (aname in anames) {
    asi <- assay_sample_info(ssamples, aname)
    
    cname <- paste0(prefix, aname)
    expect_lt(nrow(asi), nrow(has))
    sids <- has$sample_id[has[[cname]]]
    expect_setequal(sids, asi$sample_id)
  }
})

test_that("filter_by_assay_support drops samples with no assay support", {
  asi <- assay_sample_info(ssamples, "scrnaseq", drop_samples = FALSE)
  has.sc <- !is.na(asi$assay) & asi$assay == "scrnaseq"
  no.sc <- !has.sc
  
  scsamples <- filter_by_assay_support(ssamples, "scrnaseq")
  expect_setequal(scsamples$sample_id, asi$sample_id[has.sc])
  
  dropped <- samples(scsamples, dropped = TRUE)
  expect_setequal(dropped$sample_id, asi$sample_id[!has.sc])
})
