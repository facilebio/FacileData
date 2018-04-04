library(survival)
library(testthat)

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
  expect_equal(x$colnames, expected$colnames, info = varname)
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
      colnames=c(time="tte_OS", event="event_OS"),
      class="right_censored",
      label="Overall survival",

      type="clinical",
      description="Overall Survival in months"
    ))

  lol <- eav_metadata_create(pdat, covariate_def = covdef)
  fn <- tempfile()
  yaml::write_yaml(lol, fn)
  relol <- yaml::read_yaml(fn)

  # Explicitly test that the tte_OS and event_OS columns from `pDat` were
  # compounded into the OS covariatel.
  # Reference the "Encoding Survival Covariates" section in the
  # `?eav_metadata_create` helpf file for what the expected behavior of how this
  # compounded, multi-column-to-single-value mapping should work.
  compounded <- c("tte_OS", "event_OS")
  expect_true(all(compounded %in% names(pdat))) # in pData
  expect_true(!any(c("tte_OS", "event_OS") %in% names(relol))) # not in yaml
  expect_true(setequal(relol$OS$colnames, compounded)) # names of columns saved for posterity

  # ensure that variables from encoded yaml file match test meta.yaml file
  expect_true(setequal(names(relol), names(lol)))

  # encodings match
  for (varname in names(lol)) {
    validate_eav_recode(relol[[varname]], elol[[varname]], varname)
  }

  # pData with Surv
  df = data.frame(
      dataset = "foo",
      sample_id = letters[1:3],
      x = Surv(1:3, c(0,1,0)),
      y = 4:6,
      stringsAsFactors = FALSE
  )
  long = as.EAVtable(df)
  long2 = data.frame(
      dataset = "foo",
      sample_id = c("a","b","c","a","b","c"),
      variable = c("x","x","x","y","y","y"),
      value = c("1+","2 ","3+","4","5","6"),
      class = c("Surv","Surv","Surv","real","real","real"),
      stringsAsFactors = FALSE
  )
  expect_identical(long,long2)

})

test_that("basic encoding and decoding of EAV columns works", {
    # survival::Surv
    x = Surv(1:3, c(0,1,0))
    y = eav_encode_Surv(x)
    y1 = c("1+","2 ","3+")
    attr(y1, "eavclass") = "Surv"
    expect_identical(y, y1)
    z = eav_decode_Surv(y)
    expect_identical(x,z)
})
