context("Fetching assay level data")

FDS <- exampleFacileDataSet()
samples <- sample_covariate_tbl(FDS) %>%
  filter(variable == 'stage' & value == 'III') %>%
  select(dataset, sample_id)

genes <- c(
  PRF1='5551',
  GZMA='3001',
  CD274='29126',
  TIGIT='201633')

features <- tibble(assay='rnaseq', feature_id=genes)

test_that("fetch_assay_data(..., assay='rnaseq') == fetch_expression(...)", {
  ## Raw Data
  dat <- FDS %>%
    fetch_assay_data(features, samples, normalized=FALSE) %>%
    arrange(sample_id, feature_name) %>%
    select(dataset, sample_id, feature_id, symbol=feature_name, value)

  cnt <- FDS %>%
    fetch_expression(samples, feature_ids=features$feature_id) %>%
    arrange(sample_id, symbol) %>%
    select(dataset, sample_id, feature_id, symbol, value=count)

  expect_equal(dat, cnt)

  ## Normalized
  ndat <- FDS %>%
    fetch_assay_data(features, samples, normalized=TRUE) %>%
    arrange(sample_id, feature_name) %>%
    select(dataset, sample_id, feature_id, symbol=feature_name, value)

  ncnt <- FDS %>%
    fetch_expression(samples, feature_ids=features$feature_id) %>%
    cpm(prior.count=5, log=TRUE) %>%
    arrange(sample_id, symbol) %>%
    select(dataset, sample_id, feature_id, symbol, value=cpm)

  expect_equal(ndat$value, ncnt$value, tolerance=0.05)
})

test_that("fetch_assay_data(..., aggregate.by='ewm') provides scores", {
  scores <- FDS %>%
    fetch_assay_data(features, samples, normalized=TRUE, aggregate.by='ewm') %>%
    arrange(sample_id, feature_name) %>%
    select(dataset, sample_id, feature_id, symbol=feature_name, value) %>%
    mutate(samid=paste(dataset, sample_id, sep="__"))

  library(multiGSEA)
  dat <- FDS %>%
    fetch_assay_data(features, samples, normalized=TRUE, as.matrix=TRUE)
  ewm <- eigenWeightedMean(dat)$score[scores$samid]
  expect_equal(scores$value, unname(ewm))
})

test_that("fetch_assay_data handles missing entries for requested samples", {
  ## When we have multiple assays for an FDS, we can use a valid sample
  ## descriptor to retrieve data, but the requested assay may not have data
  ## for all requested samples, we need to handle this.
  root <- rprojroot::find_root(rprojroot::is_r_package)
  devtools::load_all(root)
  tcga <- FacileDataSet('~/workspace/data/facile/FacileDataSets/FacileTCGADataSet-2017-03-25')

  library(reshape2)
  samples <- sample_info_tbl(tcga) %>%
    filter(dataset == 'BRCA') %>%
    collect

  genes <- c(TIGIT='201633', CD274='29126')
  rnaseq <- tcga %>%
    fetch_assay_data(genes, samples, 'rnaseq', normalized=TRUE)

  ## don't have agilent data for all brca samples
  agilent <- tcga %>%
    fetch_assay_data(genes, samples, 'agilent', normalized=TRUE)

})
