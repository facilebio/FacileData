context("Fetching assay level data")

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
    fetch_assay_data(genes, samples, normalized = TRUE) %>%
    select(dataset, sample_id, feature_name, value) %>%
    tidyr::spread(feature_name, value)
  result <- samples %>%
    with_assay_data(genes, normalized = TRUE, .fds = FDS) %>%
    collect

  expect_equal(result, expected)
})

test_that("fetch_assay_data(..., aggregate = TRUE) provides scores", {
  scores <- FDS %>%
    fetch_assay_data(features, samples, normalized = TRUE, aggregate = TRUE) %>%
    arrange(sample_id, feature_name) %>%
    select(dataset, sample_id, feature_id, symbol=feature_name, value) %>%
    mutate(samid=paste(dataset, sample_id, sep="__"))

  dat <- FDS %>%
    fetch_assay_data(features, samples, normalized = TRUE, as.matrix = TRUE)
  ewm <- multiGSEA::eigenWeightedMean(dat)$score[scores$samid]
  expect_equal(scores$value, unname(ewm))

  # test with_assay_data
  with.scores <- scores %>%
    distinct(dataset, sample_id) %>%
    with_assay_data(features, aggregate.by = "ewm")

  expect_equal(with.scores$score,)
})

# Batch Effect Correction ======================================================

test_that("batch effect correction mimics limma::removeBatchEffect", {
  smpls <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    with_sample_covariates() %>%
    mutate(sample_key = paste(dataset, sample_id, sep = "__"))

  dat <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                          as.matrix = TRUE)
  expect_equal(smpls$sample_key, colnames(dat))

  # Normalize by one batch covaraite ...........................................
  dat.norm1 <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                                as.matrix = TRUE, batch = "sex",
                                maintain.rowmeans = FALSE)
  e.norm1 <- limma::removeBatchEffect(dat, batch = smpls$sex)
  expect_equal(dat.norm1, e.norm1)

  # the fixed expression matrix has shifted rowMeans. ..........................
  expect_true(!isTRUE(all.equal(rowMeans(dat.norm1), rowMeans(dat))))

  # We can maintain the rowmeans by setting saintain.rowmanes = TRUE
  # which is the default.
  same.mean <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                                as.matrix = TRUE, batch = "sex",
                                maintain.rowmeans = TRUE)
  expect_equal(rowMeans(same.mean), rowMeans(dat))

  # Normalize by two batch covaraites ..........................................
  set.seed(123)
  smpls$dummy <- sample(c("a", "b"), nrow(smpls), replace = TRUE)
  dat.norm2 <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                                as.matrix = TRUE, batch = c("sex", "dummy"),
                                maintain.rowmeans = FALSE)
  e.norm2 <- limma::removeBatchEffect(dat, batch = smpls$sex,
                                      batch2 = smpls$dummy)
  expect_equal(dat.norm2, e.norm2)

  # Normalize with a real valued covariate .....................................
  smpls$real <- rnorm(nrow(smpls), mean = 0)
  smpls$real[1:7] <- rnorm(7, mean = 1)
  dat.normR <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                                as.matrix = TRUE, batch = c("real"),
                                maintain.rowmeans = FALSE)
  e.normR <- limma::removeBatchEffect(dat, covariates = smpls$real)
  expect_equal(dat.normR, e.normR)

  # full bore ..................................................................
  dat.norm.uber <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                                    as.matrix = TRUE, batch = c("sex", "real"),
                                    main = "sample_type",
                                    maintain.rowmeans = FALSE)

  des <- model.matrix(~ sample_type + sex + real, smpls)
  e.norm.uber <- limma::removeBatchEffect(dat, design = des[, 1:2],
                                          covariates = des[, -(1:2)])
  expect_equal(dat.norm.uber, e.norm.uber)

  # Final check for equal rowmeans functionality ...............................
  expect_true(!isTRUE(all.equal(rowMeans(dat.norm.uber), rowMeans(dat))))
  d2 <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                         as.matrix = TRUE, batch = c("sex", "real"),
                         main = "sample_type",
                         maintain.rowmeans = TRUE)
  expect_equal(rowMeans(d2), rowMeans(dat))
})

test_that("single-gene batch correction is equivalent to all data correction", {
  set.seed(0xBEEF)
  smpls <- filter_samples(FDS, indication == "BLCA") %>%
    collect() %>%
    mutate(real.batch = rnorm(nrow(.)))

  bc.all <- fetch_assay_data(smpls, normalized = TRUE, as.matrix = TRUE,
                             batch = c("sex", "real.batch"))
  fname <- rownames(bc.all)[3]
  bc.1 <- fetch_assay_data(smpls, features = fname, normalized = TRUE,
                           as.matrix = TRUE, batch = c("sex", "real.batch"))
  expect_equal(bc.1[1,], bc.all[fname,])

  # ensure we were testing a batch corrected thing
  orig.1 <- fetch_assay_data(smpls, features = fname, normalized = TRUE,
                             as.matrix = TRUE)
  expect_equal(colnames(bc.1), colnames(orig.1))
  expect_true(!isTRUE(all.equal(bc.1[1,], orig.1[1,])))
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
