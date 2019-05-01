context("Feature Info")

FDS <- exampleFacileDataSet()
genes <- tibble(
  feature_id = c("800", "1009", "1289", "50509", "2191", "2335", "5159"),
  feature_type = "entrez")

test_that("with_feature_info grabs the right goods", {
  finfo.all <- collect(fetch_feature_info(FDS, "entrez"), n = Inf)
  finfo <- with_feature_info(genes, .fds = FDS)

  expected <- left_join(genes, finfo.all, by = c("feature_id", "feature_type"))
  expect_equal(finfo, expected)

  f2 <- with_feature_info(genes, c("name", "meta"), .fds = FDS)
  expect_equal(f2, select(expected, !!colnames(f2)))
})

test_that("with_feature_info can rename feature covariates", {
  expected <- genes %>%
    with_feature_info(c("name", "meta"), .fds = FDS) %>%
    rename(symbol = "name")

  res <- genes %>%
    with_feature_info(c(symbol = "name", "meta"), .fds = FDS)

  expect_equal(res, expected)
})
