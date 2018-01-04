context("Entity-Attribute-Value conversions")

# Checks the yaml encoding for the variable is as expected by using the
# yaml encoding from "testdata/expected-meta.yaml" matches the encoding
# that was programmatically generated
#
# @param x the recoded covariate list
# @param expected the covariate list from `expected-meta.yaml`
validate_eav_recode <- function(x, expected, varname) {
  expect_is(x, "list", info = varname)
  expect_is(expected, "list", info = varname)
  expect_equal(x$varname, expected$varname, info = varname)
  expect_equal(x$class, expected$class, info = varname)
  if (!is.null(expected$levels)) {
    expect_true(is.character(x$levels), info = varname)
    expect_equal(x$levels, expected$levels, info = varname)
  } else {
    expect_true(is.null(x$levels), info = varname)
  }
}

test_that("simple pData -> yaml covariate encoding works (no OS)", {
  # Trying to recode the survival stuff isn't included in this test
  pdat <- example_sample_covariates()
  elol <- example_sample_covariate_definitions()

  # remove OS covariates from pData and yaml
  pdat <- select(pdat, -tte_OS, -event_OS)
  elol$OS <- NULL

  lol <- create_eav_metadata(pdat)
  fn <- tempfile()
  yaml::write_yaml(lol, fn)
  relol <- yaml::read_yaml(fn)

  # ensure that all variables from dataframe were by covdefs_from_df
  expected.vars <- setdiff(colnames(pdat), c("dataset", "sample_id"))
  expect_true(setequal(names(relol), expected.vars))
  expect_true(setequal(names(relol), names(elol)))

  # encodings match
  for (varname in expected.vars) {
    validate_eav_recode(relol[[varname]], elol[[varname]], varname)
  }
})

test_that("compound pData -> meta.yaml encoding (with OS) works", {
  # Trying to recode the survival stuff isn't included in this test
  pdat <- example_sample_covariates()
  pdat <- select(pdat, sex, age, tte_OS, event_OS)
  elol <- example_sample_covariate_definitions()
  elol <- elol[c('sex','age', 'OS')]

  covdef <- list(
    OS=list(
      varname=c("tte_OS", "event_OS"),
      class="right_censored",
      label="Overall survival",

      type="clinical",
      description="Overall Survival in months"
    ))

  lol <- create_eav_metadata(pdat, covariate_def = covdef)
  fn <- tempfile()
  yaml::write_yaml(lol, fn)
  relol <- yaml::read_yaml(fn)

})
test_that("single ExpressionSet converts to FacileDataSet", {

})

test_that("list of ExpressionSets convert to FacileDataSet", {

})
