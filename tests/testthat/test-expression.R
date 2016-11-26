context("Expression")

FDS <- exampleFacileDataSet()
samples <- sample_covariate_tbl(FDS) %>%
  filter(value == 'CMS4') %>%
  select(dataset, sample_id)
samples <- sample_covariate_tbl(FDS) %>%
  filter(variable == 'stage' & value == 'III') %>%
  select(dataset, sample_id)
genes <- local({
  out <- c("800", "1009", "1289", "50509", "2191", "2335", "5159")
  intersect(out, collect(gene_info_tbl(FDS))$feature_id)
})


test_that("fetch_expression limits samples correctly", {
  s.df <- collect(samples, n=Inf)

  e.sqlite <- fetch_expression(FDS, samples, genes) %>% collect(n=Inf)
  e.df <- fetch_expression(FDS, s.df, genes) %>% collect(n=Inf)

  ## results are same from tbl_df and tbl_sqlite `samples` parameter
  expect_equal(e.sqlite, e.df)

  ## samples limited correcly
  expect_true(setequal(paste0(e.df$dataset, e.df$sample_id),
                       paste0(s.df$dataset, s.df$sample_id)))
})

test_that("fetch_expression results converted to DGEList", {
  e <- fetch_expression(FDS, samples, genes)
  y <- as.DGEList(e)
  expect_is(y, 'DGEList')

  ## check samples
  expect_is(y$samples, 'data.frame')
  expect_true(setequal(y$samples$sample_id, collect(samples)$sample_id))
  expect_type(y$samples$norm.factors, 'double')
  expect_type(y$samples$lib.size, 'integer')
  expect_type(y$samples$dataset, 'character')
  expect_type(y$samples$sample_id, 'character')

  expect_is(y$genes, 'data.frame')
  expect_type(y$genes$feature_id, 'character')
  expect_true(setequal(y$genes$feature_id, genes))
  expect_type(y$genes$symbol, 'character')
})

test_that('as.DGEList assigns correct covariates', {
  e <- fetch_expression(FDS, samples, genes)
  y <- as.DGEList(e, covariates=c('stage', 'sex'))

  expect_is(y, 'DGEList')
  expect_is(y$samples$sex, 'factor')
  expect_is(y$samples$stage, 'character')
})

test_that("cpm on fetch_expression result mimics cpm.DGEList", {
  e <- fetch_expression(FDS, samples, genes)

  E <- cpm(e, log=FALSE, prior.count=5) %>%
    mutate(samid=paste(dataset, sample_id, sep='_')) %>%
    acast(feature_id ~ samid, value.var='cpm')

  EL <- cpm(e, log=TRUE, prior.count=5) %>%
    mutate(samid=paste(dataset, sample_id, sep='_')) %>%
    acast(feature_id ~ samid, value.var='cpm')

  Etbldf <- collect(e, n=Inf) %>%
    cpm(log=TRUE, prior.count=5, .fds=FDS) %>%
    mutate(samid=paste(dataset, sample_id, sep='_')) %>%
    acast(feature_id ~ samid, value.var='cpm')

  ## counts from within database or collect should match
  expect_equal(EL, Etbldf)

  ## Test that cpm's returned with edgeR only code match
  library(edgeR)
  y <- as.DGEList(e)
  YL <- cpm(y, log=TRUE, prior.count=5)
  Y <- cpm(y, log=FALSE, prior.count=5)

  E <- E[rownames(Y), colnames(Y)]
  EL <- EL[rownames(YL), colnames(YL)]

  expect_equal(as.matrix(E), as.matrix(Y))
  expect_equal(EL, YL)
})

test_that("rpkm on fetch_expression result mimics rpkm.DGEList", {
  e <- fetch_expression(FDS, samples, genes)
  E <- rpkm(e, log=FALSE, prior.count=5) %>%
    mutate(samid=paste(dataset, sample_id, sep='_')) %>%
    acast(feature_id ~ samid, value.var='rpkm')

  EL <- rpkm(e, log=TRUE, prior.count=5) %>%
    mutate(samid=paste(dataset, sample_id, sep='_')) %>%
    acast(feature_id ~ samid, value.var='rpkm')

  Etbldf <- collect(e, n=Inf) %>%
    rpkm(log=TRUE, prior.count=5, .fds=FDS) %>%
    mutate(samid=paste(dataset, sample_id, sep='_')) %>%
    acast(feature_id ~ samid, value.var='rpkm')

  ## counts from within database or collect should match
  expect_equal(EL, Etbldf)

  ## Test that cpm's returned with edgeR only code match
  y <- as.DGEList(e)
  library(edgeR)
  Y <- rpkm(y, gene.length=y$genes$length, log=FALSE, prior.count=5)
  YL <- rpkm(y, gene.length=y$genes$length, log=TRUE, prior.count=5)

  E <- E[rownames(Y), colnames(Y)]
  EL <- EL[rownames(YL), colnames(YL)]

  expect_equal(E, Y)
  expect_equal(EL, YL)
})

