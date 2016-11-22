context("Basic FacileDataSet functions")

test_that("FacileData Constructor functions on well structure directory", {
  fds <- exampleFacileDataSet()
  genes <- c(CD96='10225', TIGIT='201633')
  exprs <- fetch_expression(fds, feature_ids=genes)
  expect_true(is(exprs, 'FacileExpression'))
})

test_that("Fetching various database tables from FacileDataSet", {
  fds <- exampleFacileDataSet()

  etable <- expression_tbl(fds)
  expect_true(is(etable, 'tbl'))

  sctable <- sample_covariate_tbl(fds)
  expect_true(is(sctable, 'tbl'))

  sstable <- sample_stats_tbl(fds)
  expect_true(is(sstable, 'tbl'))

  gitable <- gene_info_tbl(fds)
  expect_true(is(gitable, 'tbl'))
})
