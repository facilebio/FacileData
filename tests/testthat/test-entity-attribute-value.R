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

test_that("pData -> meta.yaml covariate encoding works (simple & compound)", {
  # Trying to recode the survival stuff isn't included in this test
  pdat <- example_sample_covariates()
  elol <- example_sample_covariate_definitions()

  # define covariate_def(-inition) for the compound OS facile covariate:
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

  # Explicitly test that the tte_OS and event_OS columns from `pDat` were
  # compounded into the OS covariatel.
  # Reference the "Encoding Survival Covariates" section in the
  # `?create_eav_metadata` helpf file for what the expected behavior of how this
  # compounded, multi-column-to-single-value mapping should work.
  compounded <- c("tte_OS", "event_OS")
  expect_true(all(compounded %in% names(pdat))) # in pData
  expect_true(!any(c("tte_OS", "event_OS") %in% names(relol))) # not in yaml
  expect_true(setequal(relol$OS$varname, compounded)) # names of columns saved for posterity

  # ensure that variables from encoded yaml file match test meta.yaml file
  expect_true(setequal(names(relol), names(lol)))

  # encodings match
  for (varname in names(lol)) {
    validate_eav_recode(relol[[varname]], elol[[varname]], varname)
  }
})

