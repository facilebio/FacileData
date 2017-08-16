context("Sample Covariates")

FDS <- exampleFacileDataSet()
COV.DEF <- covariate_definitions(FDS)
samples <- sample_covariate_tbl(FDS) %>%
  filter(variable == 'stage' & value == 'III') %>%
  select(dataset, sample_id)
genes <- c("800", "1009", "1289", "50509", "2191", "2335", "5159")

test_that("fetch_sample_covariates retrieves all covariates if not specified", {
  covs <- fetch_sample_covariates(FDS, samples, custom_key='lianogls') %>%
    collect(n=Inf) %>%
    arrange(dataset, sample_id, variable, value)
  ecovs <- samples %>%
    inner_join(sample_covariate_tbl(FDS), by=c('dataset', 'sample_id')) %>%
    collect(n=Inf) %>%
    arrange(dataset, sample_id, variable, value)
  expect_equal(covs, ecovs)
})

test_that("fetch_sample_covariates::samples arg limits samples correctly", {
  vars <- c('stage', 'sex', 'OS')
  covs <- fetch_sample_covariates(FDS, samples, vars) %>% collect(n=Inf)
  s.df <- collect(samples, n=Inf)
  expected.samples <- with(s.df, paste(dataset, sample_id, sep='_')) %>% unique
  returned.samples <- with(covs, paste(dataset, sample_id, sep='_')) %>% unique
  expect_true(all(returned.samples %in% expected.samples))
})

test_that("cast_covariate converts simple variables to correct type", {
  ## Simple covariates are converted to an atomic vector type
  simple.vars <- c('sex'='factor', 'stage'='character', 'sample_type'='factor')
  for (name in names(simple.vars)) {
    eclass <- simple.vars[[name]]
    info <- sprintf('%s (%s)', name, eclass)
    vals <- fetch_sample_covariates(FDS, samples, name) %>% collect(n=Inf)
    casted <- cast_covariate(name, vals$value, .fds=FDS)
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
    vals <- fetch_sample_covariates(FDS, samples, name) %>% collect(n=Inf)
    casted <- cast_covariate(name, vals$value, .fds=FDS)
    expect_is(casted, 'data.frame', info=name)
    expect_equal(names(casted), ex.names, info=name)
  }
})


test_that("spread_covariates casts simple covariates to correct class", {
  vars <- c('sex'='factor', 'stage'='character', 'sample_type'='factor')
  wide <- FDS %>%
    fetch_sample_covariates(samples, names(vars)) %>%
    spread_covariates

  ## Test presence of columns and class converted correctly
  for (cov in names(vars)) {
    expect_is(wide[[cov]], vars[cov], info=paste("Covariate:", cov))
  }

  ## Ensure that all samples asked for are in wide result
  snames <- with(collect(samples, n=Inf), paste0(dataset, '__', sample_id))
  expect_true(setequal(rownames(wide), snames))
})

test_that("spread_covariates works with both simple and complex types", {
  simple <- c('sex', 'stage')
  complex <- c('OS', 'PFS')

  ## Check values retrieved in uber result with individual results from simple
  ## and cmoplex covariates separately
  mixed <- FDS %>%
    fetch_sample_covariates(samples, c(simple, complex)) %>%
    spread_covariates

  sc <- fetch_sample_covariates(FDS, samples, simple) %>%
    spread_covariates
  ## check values in same (named) columns are the same
  cols <- intersect(names(mixed), names(sc))
  expect_equal(mixed[, cols], sc[, cols], check.attributes=FALSE)

  cc <- fetch_sample_covariates(FDS, samples, complex) %>%
    spread_covariates
  ## Some samples don't have complex covariates -- we need to fix this so that
  ## we get NA's for samples that don't have them from fetch_sample_covariates

  mcheck <- semi_join(mixed, cc, by=c('dataset', 'sample_id'))
  ## check values in same (named) columns are the same
  cols <- intersect(names(mixed), names(cc))

  mcheck %<>% arrange(dataset, sample_id)
  cc %<>% arrange(dataset, sample_id)
  expect_equal(mcheck[, cols], cc[, cols], check.attributes=FALSE)
})


test_that('with_sample_covariates returns long input with wide covariates', {
  covs <- c('stage', 'sex')

  ## Setup the individual expression and covariate table to merge into the
  ## expected result
  exprs <- fetch_assay_data(FDS, genes, samples) %>%
    collect(n=Inf)
  wcovs <- fetch_sample_covariates(FDS, samples, covs) %>%
    spread_covariates

  expected <- exprs %>%
    left_join(wcovs, by=c('dataset', 'sample_id')) %>%
    arrange(dataset, sample_id, feature_id)

  ecovs <- fetch_assay_data(FDS, genes, samples) %>%
    with_sample_covariates(covs) %>%
    arrange(dataset, sample_id, feature_id)

  expect_equal(ecovs, expected)
})
