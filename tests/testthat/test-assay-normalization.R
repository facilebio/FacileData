context("Normalizaiton of assay data")

if (!exists("FDS")) FDS <- exampleFacileDataSet()

samples <- FDS %>%
  filter_samples(stage == "III") %>%
  select(dataset, sample_id)

genes <- c(
  PRF1='5551',
  GZMA='3001',
  CD274='29126',
  TIGIT='201633')

features <- tibble(assay='rnaseq', feature_id=genes)

test_that("Normalization of rnaseq data is equivalent to edgeR::cpm", {
  y <- edgeR::calcNormFactors(as.DGEList(samples))
  cpms <- edgeR::cpm(y, log = TRUE, prior.count = 0.25)[genes,]

  # use the lib.size and norm.factors from this subset of data
  samples. <- samples %>%
    left_join(select(y$samples, sample_id, lib.size, norm.factors),
              by = "sample_id")

  normed <- fetch_assay_data(samples., genes, normalized = TRUE,
                             prior.count = 0.25, as.matrix = TRUE)
  normed <- normed[rownames(cpms), colnames(cpms)]
  expect_equal(normed, cpms)
})

