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
