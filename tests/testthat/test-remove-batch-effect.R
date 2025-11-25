context("batch effect removal")

# Currently the remove_batch_effect functionality is only tested through its
# delegation via the fetch_assay_data(..., normalized = TRUE, batch = "xxx")
# modality.

if (!exists("FDS")) FDS <- exampleFacileDataSet()

samples <- FDS |>
  filter_samples(stage == "III") |>
  select(dataset, sample_id)

genes <- c(
  PRF1='5551',
  GZMA='3001',
  CD274='29126',
  TIGIT='201633'
)

features <- tibble(assay='rnaseq', feature_id=genes)

test_that("remove_batch_effect works like limma::removeBatchEffect", {
  smpls <- FDS |>
    filter_samples(indication == "BLCA") |>
    with_sample_covariates() |>
    mutate(sample_key = paste(dataset, sample_id, sep = "__"), .before = 1L)
  
  dat <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                          as.matrix = TRUE)
  attr(dat, "fds") <- NULL
  attr(dat, "samples_dropped") <- NULL
  attr(dat, "samples") <- NULL
  
  expect_equal(smpls$sample_key, colnames(dat))
  
  # one factor covariate normalization .........................................
  norm.sex <- remove_batch_effect(dat, smpls, batch = "sex")
  norm.sexf <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                                batch = "sex", as.matrix = TRUE)
  limma.sex <- expect_message({
    limma::removeBatchEffect(dat, batch = smpls$sex)
  }, ".*not specified")
  
  expect_equal(norm.sex, limma.sex)
  expect_equal(norm.sexf, limma.sex, check.attributes = FALSE)

  # two factor coavariate normaliation .........................................
  # note that facile handles NA values in covariates and assigns them to an
  # outgroup (for better or for worse), and limma does not do that, so let's
  # manualy fix (subtype_tcga has NA values)
  set.seed(0xBEEF)
  smpls <- mutate(smpls, b2 = sample(c("a", "b"), nrow(smpls), replace = TRUE))
  norm.2 <- remove_batch_effect(dat, smpls, batch = c("sex", "b2"))
  norm.2f <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                              batch = c("sex", "b2"), as.matrix = TRUE)
  limm.2 <- expect_message({
    limma::removeBatchEffect(
      dat,
      batch = smpls$sex,
      batch2 = smpls$b2
    )    
  }, ".*not specified")

  expect_equal(norm.2, limm.2)
  expect_equal(norm.2f, limm.2, check.attributes = FALSE)

  # include a main effect to maintain ..........................................
  fbatch <- remove_batch_effect(dat, smpls, batch = "sex", main = "sample_type")
  f2 <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                         batch = "sex", main = "sample_type", as.matrix = TRUE)
  lbatch <- limma::removeBatchEffect(
    dat,
    batch = smpls$sex,
    group = smpls$sample_type
  )
  expect_equal(fbatch, lbatch)
  expect_equal(fbatch, f2, check.attributes = FALSE)

  # use a numeric covariate ....................................................
  set.seed(0xBEEF)
  smpls <- mutate(smpls, rn = rnorm(nrow(smpls)))
  fbatch <- remove_batch_effect(dat, smpls, batch = "rn")
  fb <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                         batch = "rn", as.matrix = TRUE)
  lbatch <- expect_message({
    limma::removeBatchEffect(dat, covariates = smpls$rn)
  }, ".*not specified")
  
  expect_equal(fbatch, lbatch)
  expect_equal(fbatch, fb, check.attributes = FALSE)

  # full bore ..................................................................
  dat.norm.uber <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                                    as.matrix = TRUE, batch = c("sex", "rn"),
                                    main = "sample_type")
  e.norm.uber <- limma::removeBatchEffect(dat, group = smpls$sample_type,
                                          batch = smpls$sex,
                                          covariates = smpls$rn)
  expect_equal(dat.norm.uber, e.norm.uber, check.attributes = FALSE)
})


test_that("batch correction using 'facile' covariate works", {
  smpls <- FDS |>
    filter_samples(indication == "BLCA") |>
    with_sample_covariates() |>
    mutate(sample_key = paste(dataset, sample_id, sep = "__"))

  # Unnormalized
  unnorm <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                             as.matrix = TRUE)

  # Normalize by 'sex' when 'sex' is in the sample frame
  dat.norm1 <- fetch_assay_data(FDS, features, smpls, normalized = TRUE,
                                as.matrix = TRUE, batch = "sex")
  
  # Check that all attributes (dimensions, colnames, rownames) of corrected and
  # uncorrected data are the same ...
  expect_matrix(dat.norm1, nrows = nrow(unnorm), ncols = ncol(unnorm))
  expect_equal(rownames(dat.norm1), rownames(unnorm))
  expect_equal(colnames(dat.norm1), colnames(unnorm))
  # ... but the actual values should be different
  expect_string(all.equal(unnorm, dat.norm1), pattern = "relative difference")
  
  # Normalize by 'sex' when it's not in the sample frame
  dat.norm2 <- fetch_assay_data(FDS, features, select(smpls, -sex),
                                normalized = TRUE, as.matrix = TRUE,
                                batch = "sex")
  expect_true(all.equal(dat.norm1, dat.norm2))
})

test_that("batch effect correction will handle missing covariate levels", {
  set.seed(122)
  s <- samples |>
    mutate(bcov = sample(c("a", "b", NA), nrow(samples), replace = TRUE))
  raw <- fetch_assay_data(s, genes, normalize = TRUE, as.matrix = TRUE) 
  bc <- fetch_assay_data(s, genes, batch = "bcov", normalize = TRUE,
                         as.matrix = TRUE)
  
  # this should not be equal
  expect_string(all.equal(raw, bc), pattern = "relative difference")
  
  # manually setting NA batch variables to an outgroup should match `bc`
  s <- mutate(s, bcov2 = ifelse(is.na(bcov), "outgroup", bcov))
  bc.man <- remove_batch_effect(raw, s, batch = "bcov2")
  expect_equal(bc.man, bc)
})

test_that("single-gene batch correction is equivalent to all data correction", {
  set.seed(0xBEEF)
  smpls <- filter_samples(FDS, indication == "BLCA") |>
    collect(n = INF) |>
    mutate(real.batch = rnorm(n()))

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
