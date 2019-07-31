context("Basic FacileDataSet functions")

FDS <- exampleFacileDataSet()

test_that("Fetching various database tables from FacileDataSet", {
  sctable <- sample_covariate_tbl(FDS)
  expect_true(is(sctable, 'tbl'))

  sstable <- sample_stats_tbl(FDS)
  expect_true(is(sstable, 'tbl'))

  gitable <- gene_info_tbl(FDS)
  expect_true(is(gitable, 'tbl'))
})

test_that("compound filter criteria == method chaining with filter_samples()", {
  s1 <- FDS %>%
    filter_samples(indication == "CRC", sex == "f")
  s2 <- FDS %>%
    filter_samples(indication == "CRC") %>%
    filter_samples(sex == "f")
  assert_set_equal(s2$sample_id, s1$sample_id)
})

test_that("filter_samples filters against dataset and sample_id columns", {
  all.samples <- FDS %>%
    samples() %>%
    with_sample_covariates()

  blca.f <- filter_samples(FDS, dataset == "BLCA", sex == "f")
  blca.e <- filter(all.samples, dataset == "BLCA", sex == "f")
  assert_set_equal(blca.f$sample_id, blca.e$sample_id)

  some.ids <- sample(all.samples$sample_id, 5)
  some.f <- filter_samples(FDS, sample_id %in% some.ids)
  some.e <- filter(all.samples, sample_id %in% some.ids)
  assert_set_equal(some.f$sample_id, some.e$sample_id)
})
