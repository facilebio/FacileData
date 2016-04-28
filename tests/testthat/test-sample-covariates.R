context("Sample Covariate")

DB <- FacileDb()
samples <- sample_covariate_tbl(DB) %>%
  filter(value == 'TNBC') %>%
  select(dataset, sample_id)
# samples <- sample_covariate_tbl(DB) %>%
#   filter(value == 'IC3') %>%
#   select(dataset, sample_id)

test_that("spread_covariates makes wide data.frame", {
  s.df <- collect(samples)
  vars <- c('BCOR', 'IC', 'TC', 'OS', 'subtype_receptor')
  covs <- fetch_sample_covariates(DB, samples, vars)
  wide <- spread_covariates(covs)

  expect_is(wide$BCOR, 'factor')
  expect_is(wide$IC, 'factor')
  expect_is(wide$TC, 'factor')

  expect_is(wide$OS, 'numeric')

  ## samples limited correcly
  expect_true(setequal(rownames(wide),
                       paste0(s.df$dataset, '_', s.df$sample_id)))
})
