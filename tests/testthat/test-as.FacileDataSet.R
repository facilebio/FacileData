context("as.FacileDataSet")

es <- multiGSEA::exampleExpressionSet(do.voom=FALSE)
esl <- list(first=es, second=es)
colnames(esl[['second']]) <- paste0('two_', colnames(esl[['second']]))

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
  message("=== ", varname)
}

test_that("propper yaml covariate extraction from pData", {
  pdat <- system.file("testdata", "test-sample-covariates.rds",
                      package = "FacileDataSet")
  pdat <- readRDS(pdat)
  elol <- system.file("testdata", "expected-meta.yaml",
                      package = "FacileDataSet")
  elol <- yaml::read_yaml(elol)$sample_covariate

  lol <- covdefs_from_df(pdat)
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

test_that("single ExpressionSet converts to FacileDataSet", {

})

test_that("list of ExpressionSets convert to FacileDataSet", {

})
