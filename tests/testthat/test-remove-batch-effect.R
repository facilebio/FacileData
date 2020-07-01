context("batch effect removal")

# Currently the remove_batch_effect functionality is only tested through its
# delegation via the fetch_assay_data(..., normalized = TRUE, batch = "xxx")
# modality.

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

test_that("batch correction using 'facile' covariate works", {
  smpls <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    with_sample_covariates() %>%
    mutate(sample_key = paste(dataset, sample_id, sep = "__"))

  # Unnormalized
  unnorm <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                             as.matrix = TRUE)

  # Normalize by 'sex' when 'sex' is in the sample frame
  dat.norm1 <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                                as.matrix = TRUE, batch = "sex")
  expect_matrix(dat.norm1, nrows = nrow(unnorm), ncols = ncol(unnorm))
  expect_equal(rownames(dat.norm1), rownames(unnorm))
  expect_equal(colnames(dat.norm1), colnames(unnorm))

  # Normalize by 'sex' when it's not in the sample frame
  dat.norm2 <- fetch_assay_data(FDS, features, select(smpls, -sex),
                                normalized = TRUE, as.matrix = TRUE,
                                batch = "sex")

  expect_true(!isTRUE(all.equal(dat.norm1, unnorm)))
  expect_true(all.equal(dat.norm1, dat.norm2))
})

test_that("batch effect correction will handle missing covariate levels", {
  set.seed(122)
  s <- samples %>%
    mutate(bcov = sample(c("a", "b", NA), nrow(samples), replace = TRUE))
  bc <- fetch_assay_data(s, genes, batch = "bcov", normalize = TRUE,
                         as.matrix = TRUE)
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

test_that("including `main` with `batch` works as expected", {
  # TODO: compare batch correction using airway dataset. It seems like including
  # `main = "cell"` doesn't have any effect, as it is being performed in the
  # FacileAnalysis,FacileAnalysis-RNAseq.Rmd vignette
  #
  # Shouldn't these two be different?
  # bc.1 <- remove_batch_effect(airway, batch = "cell")
  # bc.2 <- remove_batch_effect(airway, batch = "cell", main = "dex")
})
