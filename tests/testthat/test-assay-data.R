context("Fetching assay level data")

FDS <- exampleFacileDataSet()
samples <- sample_covariate_tbl(FDS) %>%
  filter(variable == 'stage' & value == 'III') %>%
  select(dataset, sample_id)

samples <- FDS %>%
  filter_samples(stage == "III") %>%
  select(dataset, sample_id)

genes <- c(
  PRF1='5551',
  GZMA='3001',
  CD274='29126',
  TIGIT='201633')

features <- tibble(assay='rnaseq', feature_id=genes)

test_that("fetch_assay_data limits samples correctly", {
  s.df <- collect(samples, n=Inf)

  e.sqlite <- fetch_assay_data(FDS, genes, samples) %>% collect(n=Inf)
  e.df <- fetch_assay_data(FDS, genes, s.df) %>% collect(n=Inf)

  ## results are same from tbl_df and tbl_sqlite `samples` parameter
  expect_equal(e.sqlite, e.df)

  ## samples limited correcly
  expect_true(setequal(paste0(e.df$dataset, e.df$sample_id),
                       paste0(s.df$dataset, s.df$sample_id)))

})

test_that("spreading data works with_assay_data", {
  expected <- FDS %>%
    fetch_assay_data(genes, samples, normalized=TRUE) %>%
    select(dataset, sample_id, feature_name, value) %>%
    tidyr::spread(feature_name, value)
  result <- samples %>%
    with_assay_data(genes, normalized = TRUE, .fds = FDS) %>%
    collect

  expect_equal(result, expected)
})

test_that("fetch_assay_data(..., aggregate.by='ewm') provides scores", {
  scores <- FDS %>%
    fetch_assay_data(features, samples, normalized=TRUE, aggregate.by='ewm') %>%
    arrange(sample_id, feature_name) %>%
    select(dataset, sample_id, feature_id, symbol=feature_name, value) %>%
    mutate(samid=paste(dataset, sample_id, sep="__"))

  dat <- FDS %>%
    fetch_assay_data(features, samples, normalized=TRUE, as.matrix=TRUE)
  ewm <- multiGSEA::eigenWeightedMean(dat)$score[scores$samid]
  expect_equal(scores$value, unname(ewm))

  xx <- scores %>%
    distinct(dataset, sample_id) %>%
    with_assay_data(features, normalized = TRUE, aggregate.by = "ewm")
})

# test_that("fetch_assay_data handles missing entries for requested samples", {
#   ## When we have multiple assays for an FDS, we can use a valid sample
#   ## descriptor to retrieve data, but the requested assay may not have data
#   ## for all requested samples, we need to handle this.
#   root <- rprojroot::find_root(rprojroot::is_r_package)
#   devtools::load_all(root)
#   tcga <- FacileDataSet('~/workspace/data/facile/FacileDataSets/FacileTCGADataSet-2017-03-25')
#
#   library(reshape2)
#   samples <- sample_info_tbl(tcga) %>%
#     filter(dataset == 'BRCA') %>%
#     collect
#
#   genes <- c(TIGIT='201633', CD274='29126')
#   rnaseq <- tcga %>%
#     fetch_assay_data(genes, samples, 'rnaseq', normalized=TRUE)
#
#   ## don't have agilent data for all brca samples
#   agilent <- tcga %>%
#     fetch_assay_data(genes, samples, 'agilent', normalized=TRUE)
#
# })
