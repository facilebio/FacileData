context("EAV Manipulation")

efds. <- exampleFacileDataSet()
pdata. <- efds. %>%
  # sex and stage are factors
  fetch_sample_covariates(covariates = c("sex", "stage")) %>%
  spread_covariates() %>%
  # add some other data types: numeric and character
  mutate(age = sample(20:70, nrow(.)),
         category = sample(letters, nrow(.), replace = TRUE))
emeta. <- local({
  fn <- system.file("extdata", "exampleFacileDataSet", "meta.yaml",
                    package = "FacileData", mustWork = TRUE)
  defined <- yaml::yaml.load_file(fn)$sample_covariates
  defined <- defined[names(defined) %in% colnames(pdata.)]
  c(defined, list(
    age = list(type = "atype", class = "real", description = "happy bday"),
    category = list(type = "atype", class = "categorical", description = "x")))
})

test_that("deafult metadata creation from data.frame", {
  ignore.cols <- c("dataset", "sample_id")

  covdefs <- eav_metadata_create(pdata., ignore = ignore.cols)
  expected.cols <- setdiff(names(pdata.), ignore.cols)
  expect_setequal(names(covdefs), expected.cols)

  # check that inferred covariate definitions have the required "slots", ie.
  # arguments, label, class, type
  for (cname in expected.cols) {
    vals <- pdata.[[cname]]
    expected <- emeta.[[cname]]
    inferred <- covdefs[[cname]]
    if (is.character(vals) || is.factor(vals)) {
      expect_equal(inferred[["class"]], "categorical", info = cname)
    } else if (is.numeric(vals)) {
      expect_equal(inferred[["class"]], "real", info = cname)
    }
    if (is.factor(vals)) {
      expect_equal(inferred[["levels"]], levels(vals), info = cname)
    }
  }
})

test_that("custom definition supersede inferred EAV defs from data.frame", {
  covariate_def <- list(
    # Change order of levels and description in sex factor covariate
    sex = list(
      levels = rev(levels(pdata.[["sex"]])),
      description = "x:y chromosome ratio"),
    age = list(
      description = "years of life",
      type = "random"))

  defaults <- eav_metadata_create(pdata.)
  custom <-  eav_metadata_create(pdata., covariate_def = covariate_def)

  # Check reversed levels of `sex` covariate
  expect_character(custom$sex$levels)
  expect_equal(custom$sex$levels, rev(defaults$sex$levels))

  # Ensure that specified entries were overriden, and others left alone.
  for (cname in names(defaults)) {
    dvals <- defaults[[cname]]
    cvals <- custom[[cname]]
    assert_list(dvals)
    assert_list(cvals)
    assert_subset(names(dvals), names(cvals))
    for (attrib in names(dvals)) {
      override <- covariate_def[[cname]][[attrib]]
      expected <- if (!is.null(override)) override else dvals[[attrib]]
      expect_equal(cvals[[attrib]], expected,
                   info = paste(cname, attrib, sep = ":"))
    }
  }
})
