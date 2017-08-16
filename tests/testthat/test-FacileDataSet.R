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
