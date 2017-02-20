context("Expression")

FDS <- exampleFacileDataSet()
samples <- sample_covariate_tbl(FDS) %>%
  filter(variable == 'stage' & value == 'III') %>%
  select(dataset, sample_id)
genes <- local({
  out <- c("800", "1009", "1289", "50509", "2191", "2335", "5159")
  feature_info_tbl(FDS) %>%
    filter(feature_id %in% out) %>%
    collect %$%
    feature_id
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

test_that("spread_expression works", {
  exprs <- fetch_expression(FDS, samples, genes) %>%
    cpm(log=TRUE, prior.count=5) %>%
    arrange(dataset, sample_id, feature_id)

  sfcnt <- spread_expression(exprs, 'feature_id', 'count')
  sfcpm <- spread_expression(exprs, 'feature_id', 'cpm')

  ## Check that colnames of spread count and cpm expression have all feature_ids
  expected_fids <- unique(exprs$feature_id)
  fids.cnt <- local({
    tmp <- colnames(sfcnt)[grep('feature_id', colnames(sfcnt))]
    sub('feature_id_', '', tmp)
  })
  fids.cpm <- local({
    tmp <- colnames(sfcpm)[grep('feature_id', colnames(sfcpm))]
    sub('feature_id_', '', tmp)
  })
  expect_true(setequal(fids.cnt, expected_fids))
  expect_true(setequal(fids.cpm, expected_fids))

  ## Check that the columns of the spread features are of the correct type
  for (fid in expected_fids) {
    cnts <- sfcnt[[paste0('feature_id_', fid)]]
    cpms <- sfcpm[[paste0('feature_id_', fid)]]
    expect_true(is(cnts, 'integer'), info=fid)
    expect_true(is(cpms, 'numeric'), info=fid)
  }

  ## Check that the spread counts match original values
  sf.counts <- sfcnt %>%
    melt(id.vars=c('dataset', 'sample_id'),
         variable.name='feature_id', value.name='count') %>%
    mutate(feature_id=sub('feature_id_', '', feature_id)) %>%
    arrange(dataset, sample_id, feature_id)
  expect_equal(sf.counts, select(exprs, dataset:count), check.attributes=FALSE)

  ## Less rigorous: check that columns of spread symbols have all the symbols
  ## in the column names that existed in the original data.
  sscnt <- spread_expression(exprs, 'symbol', 'count')
  sscpm <- spread_expression(exprs, 'symbol', 'cpm')
  expect_true(length(setdiff(exprs$symbol, colnames(sscnt))) == 0)
  expect_true(length(setdiff(exprs$symbol, colnames(sscpm))) == 0)
})

test_that("fetch_expression results converted to DGEList", {
  e <- fetch_expression(FDS, samples, genes)
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
})

test_that('as.DGEList assigns correct covariates', {
  e <- fetch_expression(FDS, samples, genes)
  y <- as.DGEList(e, covariates=c('indication', 'sex'))

  expect_is(y, 'DGEList')
  expect_is(y$samples$sex, 'factor')
  expect_is(y$samples$indication, 'character')
})

test_that("as.DGEList accepts character or covariate data.frame", {
  cov.df <- fetch_sample_covariates(FDS, samples, c('sex', 'indication'))
  y0 <- as.DGEList(samples, c('sex', 'indication'), .fds=FDS)
  y <- as.DGEList(samples, cov.df, .fds=FDS)
  expect_equal(y$samples, y0$samples)
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

