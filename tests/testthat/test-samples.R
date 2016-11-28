context("Retrieving arbitrary samples")

FDS <- exampleFacileDataSet()

test_that("fetch_samples allows filtering covariate table as if it were wide", {
  samples <- FDS %>%
    fetch_samples(sex == 'f', subtype_crc_cms %in% c('CMS1', 'CMS2'))
  asamples <- with_sample_covariates(samples, c('sex', 'indication', 'subtype_crc_cms'))
  asamples <- droplevels(asamples)

  ## This should only retrieve READ and COAD datasets
  expect_true(setequal(asamples$dataset, c("READ", "COAD")))
  expect_true(all(asamples$sex == 'f'))
  expect_true(setequal(asamples$subtype_crc_cms, c('CMS1', 'CMS2')))
})

test_that("fetch_samples throws error when subsetting with unknown covariate", {
  asamples <- FDS %>%
    fetch_samples(stage == 'I') %>%
    with_sample_covariates('stage')
  expect_true(all(asamples$stage == 'I'))

  expect_error(FDS %>% fetch_samples(stages == 'I'), 'not defined: stages')
})
