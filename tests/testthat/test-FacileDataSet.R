context("Basic FacileDataSet functions")

FDS <- exampleFacileDataSet()

test_that("FacileData Constructor functions on well structure directory", {
  genes <- c(CD96='10225', TIGIT='201633')
  exprs <- fetch_expression(FDS, feature_ids=genes)
  expect_true(is(exprs, 'FacileExpression'))
})

test_that("Fetching various database tables from FacileDataSet", {
  sctable <- sample_covariate_tbl(FDS)
  expect_true(is(sctable, 'tbl'))

  sstable <- sample_stats_tbl(FDS)
  expect_true(is(sstable, 'tbl'))

  gitable <- gene_info_tbl(FDS)
  expect_true(is(gitable, 'tbl'))
})
