context("Sample Covariates")

FDS <- exampleFacileDataSet()
COV.DEF <- covariate_definitions(FDS)
samples <- sample_covariate_tbl(FDS) |>
  filter(variable == 'stage' & value == 'III') |>
  distinct(dataset, sample_id)
genes <- c("800", "1009", "1289", "50509", "2191", "2335", "5159")

test_that("fetch_sample_covariates queries all samples if not specified", {
  covs <- fetch_sample_covariates(FDS, "sample_type")
  expect_setequal(covs[["sample_id"]], collect(samples(FDS))[["sample_id"]])
})

test_that("fetch_sample_covariates retrieves all covariates if not specified", {
  covs <- fetch_sample_covariates(FDS, samples = samples) |>
    collect(n = Inf) |>
    arrange(dataset, sample_id, variable, value)
  ecovs <- samples |>
    inner_join(sample_covariate_tbl(FDS), by = c('dataset', 'sample_id')) |>
    collect(n = Inf) |>
    arrange(dataset, sample_id, variable, value)
  expect_equal(covs, ecovs, check.attributes = FALSE)
})

test_that("fetch_sample_covariates::samples arg limits samples correctly", {
  vars <- c('stage', 'sex', 'OS')
  covs <- fetch_sample_covariates(FDS, vars, samples) |> collect(n = Inf)
  s.df <- collect(samples, n=Inf)
  expected.samples <- with(s.df, paste(dataset, sample_id, sep='_')) |> unique()
  returned.samples <- with(covs, paste(dataset, sample_id, sep='_')) |> unique()
  expect_true(all(returned.samples %in% expected.samples))
})

test_that("cast_covariate converts simple variables to correct type", {
  ## Simple covariates are converted to an atomic vector type
  simple.vars <- c(
    sex = "factor",
    stage ="factor", sample_type = "factor",
    subtype_tcga = "character",
    OS = "data.frame")

  for (name in names(simple.vars)) {
    eclass <- simple.vars[[name]]
    info <- sprintf('%s (%s)', name, eclass)
    vals <- fetch_sample_covariates(FDS, name, samples) |> collect(n = Inf)
    casted <- cast_covariate(name, vals$value, .fds = FDS)
    len <- if (eclass == "data.frame") nrow(casted) else length(casted)
    expect_true(len == nrow(vals), info = info)
    expect_is(casted, eclass, info = info)
  }
})

test_that("cast_covariate converts right_censored data correctly", {
  ## right_censored data should be converted into a 2-column data.frame with
  ## propper names
  vars <- c('OS', 'PFS')
  for (name in vars) {
    ex.names <- paste(c('tte', 'event'), name, sep='_')
    vals <- fetch_sample_covariates(FDS, name, samples) |> collect(n = Inf)
    casted <- cast_covariate(name, vals$value, .fds = FDS)
    expect_is(casted, 'data.frame', info = name)
    expect_true(setequal(names(casted), ex.names), info = name)
    expect_is(casted[[ex.names[1]]], "numeric", info = ex.names[1])
    expect_is(casted[[ex.names[2]]], "integer", info = ex.names[2])
  }
})


test_that("spread_covariates casts simple covariates to correct class", {
  vars <- c(
    sex = "factor",
    stage ="factor", sample_type = "factor",
    subtype_tcga = "character")

  wide <- FDS |>
    fetch_sample_covariates(names(vars), samples) |>
    spread_covariates()

  ## Test presence of columns and class converted correctly
  for (cov in names(vars)) {
    expect_is(wide[[cov]], vars[cov], info = paste("Covariate:", cov))
  }

  ## Ensure that all samples asked for are in wide result
  missed <- anti_join(wide, collect(samples), by = c('dataset', 'sample_id'))
  expect_true(nrow(missed) == 0)
})

test_that("spread_covariates works with both simple and complex types", {
  simple <- c('sex', 'stage')
  complex <- c('OS', 'PFS')

  # Check values retrieved in uber result with individual results from simple
  # and cmoplex covariates separately
  mixed <- FDS |>
    fetch_sample_covariates(c(simple, complex), samples) |>
    spread_covariates()

  sc <- fetch_sample_covariates(FDS, simple, samples) |>
    spread_covariates()
  # check values in same (named) columns are the same
  cols <- intersect(names(mixed), names(sc))
  expect_equal(mixed[, cols], sc[, cols], check.attributes=FALSE)

  cc <- fetch_sample_covariates(FDS, complex, samples) |>
    spread_covariates()
  # Some samples don't have complex covariates -- we need to fix this so that
  # we get NA's for samples that don't have them from fetch_sample_covariates

  mcheck <- semi_join(mixed, cc, by = c("dataset", "sample_id"))
  # check values in same (named) columns are the same
  cols <- intersect(names(mixed), names(cc))

  mcheck <- arrange(mcheck, dataset, sample_id)
  cc <- arrange(cc, dataset, sample_id)
  expect_equal(mcheck[, cols], cc[, cols], check.attributes = FALSE)
})


test_that('with_sample_covariates returns long input with wide covariates', {
  covs <- c('stage', 'sex')

  ## Setup the individual expression and covariate table to merge into the
  ## expected result
  exprs <- fetch_assay_data(FDS, genes, samples) |>
    collect(n = Inf)
  wcovs <- fetch_sample_covariates(FDS, covs, samples) |>
    spread_covariates()

  expected <- exprs |>
    left_join(wcovs, by=c('dataset', 'sample_id')) |>
    arrange(dataset, sample_id, feature_id)

  ecovs <- fetch_assay_data(FDS, genes, samples) |>
    with_sample_covariates(covs) |>
    arrange(dataset, sample_id, feature_id)

  expect_equal(ecovs, expected, check.attributes = FALSE)
})

test_that("successive with_sample_covariate calls build correct frame", {
  expected <- samples |> 
    fetch_sample_covariates(c('sex', 'stage')) |> 
    spread_covariates() |>
    arrange(dataset, sample_id)

  res <- samples |>
    with_sample_covariates("sex") |>
    with_sample_covariates("stage") |>
    arrange(dataset, sample_id)
  expect_equal(res, expected)
})

test_that("with_sample_covariates can support renaming covariates", {
  expected <- samples |>
    fetch_sample_covariates(c('sample_type', 'stage')) |>
    spread_covariates() |>
    arrange(dataset, sample_id) |>
    rename(tumor_stage = "stage")
  res <- samples |>
    with_sample_covariates(c("sample_type", tumor_stage = "stage")) |>
    arrange(dataset, sample_id)
  expect_equal(res, expected)
})

test_that("summary.[eav_covariates|wide_covariates] make congruent results", {
  FDS <- an_fds()
  asamples <- samples(FDS)
  long <- fetch_sample_covariates(asamples)
  wide <- with_sample_covariates(asamples)  
  
  lres <- summary(long, expanded = FALSE) |>
    arrange(variable)
  wres <- summary(wide, expanded = FALSE) |> 
    arrange(variable)
  expect_equal(lres, wres)  
  
  lres.exp <- summary(long, expanded = TRUE) |> 
    arrange(variable, level)
  wres.exp <- summary(wide, expanded = TRUE) |> 
    arrange(variable, level)
  expect_equal(lres.exp, wres.exp)
})

test_that("covariate summary is robust to sparse numeric vectors", {
  # when a numeric vector is short or has many NA's, the quantile & cut mojo
  # fails due to insufficient number of intervals, and what have you.
  #   
  #    Error in cut.default(vals, qtl, include.lowest = TRUE) :
  #      invalid number of intervals
  efds <- exampleFacileDataSet()
  wide <- samples(efds) |> with_sample_covariates()
  long <- samples(efds) |> fetch_sample_covariates()
  
  wpdat <- samples(efds) |>
    with_sample_covariates() |> 
    summary(expanded = TRUE) |> 
    arrange(variable, level)
  lpdat <- samples(efds) |> 
    fetch_sample_covariates() |> 
    summary(expanded = TRUE) |> 
    arrange(variable, level)
  # Note that these won't be equal because the wide format converts the 
  # "event" (response) covariates into two columns in "a smart way", but
  # the long format doesn't do jack on that.
  # 
  # This test for now just checks that no error was thrown
  expect_s3_class(wpdat, "facile_frame")
  expect_s3_class(lpdat, "facile_frame")
})
