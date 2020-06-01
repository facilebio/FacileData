context("samples(FacileDataSet)")

test_that("samples() is a facile_frame", {
  efds <- exampleFacileDataSet()

  expected <- dplyr::tbl(efds$con, 'sample_info') %>%
    collect() %>%
    select(dataset, sample_id) %>%
    arrange(sample_id)
  samples. <- samples(efds) %>%
    collect() %>%
    arrange(sample_id)
  expect_equal(samples., expected, check.attributes = FALSE)
  expect_s3_class(fds(samples.), "FacileDataSet")
})
