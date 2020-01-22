context("Testing conversion to Bioc Expression Containers")

if (!exists("FDS")) FDS <- exampleFacileDataSet()

samples <- sample_covariate_tbl(FDS) %>%
  filter(variable == 'stage' & value == 'III') %>%
  select(dataset, sample_id)
genes <- local({
  out <- c("800", "1009", "1289", "50509", "2191", "2335", "5159")
  feature_info_tbl(FDS) %>%
    filter(feature_id %in% out) %>%
    collect() %>%
    pull(feature_id)
})

test_that("fetch_assay_data results converted to DGEList", {
  e <- fetch_assay_data(FDS, genes, samples)
  y <- as.DGEList(e)
  expect_is(y, 'DGEList')

  ## check samples
  expect_is(y$samples, 'data.frame')
  expect_true(setequal(y$samples$sample_id, collect(samples)$sample_id))
  expect_type(y$samples$norm.factors, 'double')
  expect_type(y$samples$lib.size, 'double')
  expect_type(y$samples$dataset, 'character')
  expect_type(y$samples$sample_id, 'character')

  expect_is(y$genes, 'data.frame')
  expect_type(y$genes$feature_id, 'character')
  expect_true(setequal(y$genes$feature_id, genes))
  expect_type(y$genes$symbol, 'character')

  ## Check that counts match up in DGEList as they would from raw matrix fetch
  m <- fetch_assay_data(FDS, genes, samples, normalize=FALSE, as.matrix=TRUE)
  expect_true(setequal(rownames(m), rownames(y)))
  expect_true(setequal(colnames(m), colnames(y)))
  expect_equal(m[rownames(y), colnames(y)], y$counts)
})

test_that("as.DGEList appends custom covariate table correctly", {
  custom.covs <- c("sex", "subtype_molecular")
  with.covs <- with_sample_covariates(samples, custom.covs)

  y.ref <- as.DGEList(with.covs)
  y.test <- as.DGEList(with.covs, covariates = with.covs)

  expect_equal(dim(y.test), dim(y.ref))
  expect_equal(colnames(y.test), colnames(y.ref))

  expected.scols <- c(
    "group", "lib.size", "norm.factors", "dataset", "sample_id", "samid",
    custom.covs)

  expect_set_equal(colnames(y.test$samples), expected.scols)

  for (cov in custom.covs) {
    expect_equal(y.test$samples[[cov]], y.ref$samples[[cov]],
                 info = paste("sample <-> custom covariate match:", cov))
  }
})

test_that("as.DGEList with custom lib.size and norm.factors works", {
  set.seed(100)
  y.all <- as.DGEList(samples)
  y.some <- y.all[sample(nrow(y.all), 1000),, keep.lib.sizes = FALSE]
  y.some <- edgeR::calcNormFactors(y.some)

  assert_true(all(y.all$samples$lib.size > y.some$samples$lib.size))
  assert_false(any(y.all$samples$norm.factors == y.some$samples$norm.factors))

  # cpm calculation with stored libsize and normfactors
  cpms.orig <- fetch_assay_data(samples, features = rownames(y.some),
                                as.matrix = TRUE, normalized = TRUE, log = TRUE,
                                prior.count = 2)

  # Add custom lib.size and norm.factors to sample facie_frame
  samples. <- samples %>%
    collect() %>%
    left_join(select(y.some$samples, sample_id, lib.size, norm.factors),
              by = "sample_id")

  cpms.f <- fetch_assay_data(samples., features = rownames(y.some),
                             as.matrix = TRUE, normalized = TRUE, log = TRUE,
                             prior.count = 2)
  expect_equal(rownames(cpms.orig), rownames(cpms.f))
  expect_equal(colnames(cpms.orig), colnames(cpms.f))

  # This is an explicit test to make sure that the differences in the CPMs are
  # not zero
  mean.diff.orig <- mean(abs(cpms.f - cpms.orig))
  expect_gt(mean.diff.orig, 0)

  cpms.e1 <- edgeR::cpm(y.some, log = TRUE, prior.count = 2)
  expect_setequal(rownames(cpms.f), rownames(cpms.e1))
  expect_setequal(colnames(cpms.f), colnames(cpms.e1))
  cpms.f <- cpms.f[rownames(cpms.e1), colnames(cpms.e1)]
  mean.diff.new <- mean(abs(cpms.f - cpms.e1))
  expect_equal(mean.diff.new, 0)
})
