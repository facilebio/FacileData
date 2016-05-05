context("Sample Covariate")

DB <- FacileDb()
samples <- sample_covariate_tbl(DB) %>%
  filter(value == 'TNBC') %>%
  select(dataset, sample_id)
genes <- c("800", "1009", "1289", "50509", "2191", "2335", "5159")

test_that("fetch_sample_covariates::samples arg limits samples correctly", {
  vars <- c('BCOR', 'IC', 'TC', 'OS')
  covs <- fetch_sample_covariates(DB, samples, vars) %>% collect
  s.df <- collect(samples)
  expected.samples <- with(s.df, paste(dataset, sample_id, sep='_')) %>% unique
  returned.samples <- with(covs, paste(dataset, sample_id, sep='_')) %>% unique
  expect_true(all(returned.samples %in% expected.samples))
})

test_that("spread_covariates casts (with class) long result to wide & dense tbl", {
  vars <- c('BCOR', 'IC', 'TC', 'OS')
  covs <- fetch_sample_covariates(DB, samples, vars)
  wide <- spread_covariates(covs)

  ## Test presence of columns and class convertec correctly
  expect_is(wide$BCOR, 'factor')
  expect_is(wide$IC, 'factor')
  expect_is(wide$TC, 'factor')
  expect_is(wide$OS, 'numeric')

  ## Ensure that all samples asked for are in wide result
  snames <- with(collect(samples), paste0(dataset, '_', sample_id))
  expect_true(setequal(rownames(wide), snames))
})


test_that("spread_covariates(..., cov.def=NULL) sets column-values as character", {
  vars <- c('BCOR', 'IC', 'TC', 'OS', 'subtype_receptor')
  covs <- fetch_sample_covariates(DB, samples, vars)
  expected <- collect(covs) %>%
    dcast(dataset + sample_id ~ variable, value.var='value') %>%
    set_rownames(., paste0(.$dataset, '_', .$sample_id))
  result <- spread_covariates(covs, NULL)

  expect_equal(result, expected)
  expect_true(all(sapply(result, class) == 'character'))
})

test_that('with_sample_covariates returns long input with wide covariates', {
  covs <- c('BCOR', 'IC', 'TC')

  ## Setup the individual expression and covariate table to merge into the
  ## expected result
  exprs <- fetch_expression(DB, samples, genes) %>%
    collect
  wcovs <- fetch_sample_covariates(DB, samples, covs) %>%
    spread_covariates
  expected <- left_join(exprs, wcovs, by=c('dataset', 'sample_id')) %>%
    arrange(dataset, sample_id, feature_id)

  ecovs <- fetch_expression(DB, samples, genes) %>%
    with_sample_covariates(covs) %>%
    arrange(dataset, sample_id, feature_id)

  expect_equal(ecovs, expected)
})
