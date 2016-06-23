context("Sample Covariate")

DB <- TestDb()
cov.def <- DB[['cov.def']]
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

test_that("cast_covariate converts simple variables to correct type", {
  ## Simple covariates are converted to an atomic vector type
  simple.vars <- c('BCOR'='factor', 'IC'='factor', 'TC'='factor')
  for (name in names(simple.vars)) {
    eclass <- simple.vars[[name]]
    info <- sprintf('%s (%s)', name, eclass)
    vals <- fetch_sample_covariates(DB, samples, name) %>% collect
    casted <- cast_covariate(name, vals$value, cov.def)
    expect_true(length(casted) == nrow(vals), info=info)
    expect_equal(class(casted)[1L], eclass, info=info)
  }
})

test_that("cast_covariate converts right_censored data correctly", {
  ## right_censored data should be converted into a 2-column data.frame with
  ## propper names
  vars <- c('OS', 'PFS')
  for (name in vars) {
    ex.names <- paste(c('tte', 'event'), name, sep='_')
    vals <- fetch_sample_covariates(DB, samples, name) %>% collect
    casted <- cast_covariate(name, vals$value, cov.def)
    expect_is(casted, 'data.frame', info=name)
    expect_equal(names(casted), ex.names, info=name)
  }
})


test_that("spread_covariates casts simple covariates to correct class", {
  vars <- c('BCOR', 'IC', 'TC')
  wide <- fetch_sample_covariates(DB, samples, vars) %>%
    spread_covariates

  ## Test presence of columns and class convertec correctly
  expect_is(wide$BCOR, 'factor')
  expect_is(wide$IC, 'factor')
  expect_is(wide$TC, 'factor')

  ## Ensure that all samples asked for are in wide result
  snames <- with(collect(samples), paste0(dataset, '_', sample_id))
  expect_true(setequal(rownames(wide), snames))
})

test_that("spread_covariates works with both simple and complex types", {
  simple <- c('BCOR', 'IC')
  complex <- c('OS', 'PFS')

  ## Check values retrieved in uber result with individual results from simple
  ## and cmoplex covariates separately
  mixed <- fetch_sample_covariates(DB, samples, c(simple, complex)) %>%
    spread_covariates

  sc <- fetch_sample_covariates(DB, samples, simple) %>%
    spread_covariates
  ## check values in same (named) columns are the same
  cols <- intersect(names(mixed), names(sc))
  expect_equal(mixed[, cols], sc[, cols], check.attributes=FALSE)

  cc <- fetch_sample_covariates(DB, samples, complex) %>%
    spread_covariates
  ## check values in same (named) columns are the same
  cols <- intersect(names(mixed), names(cc))
  expect_equal(mixed[, cols], cc[, cols], check.attributes=FALSE)
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
