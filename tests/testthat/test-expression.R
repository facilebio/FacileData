context("Expression")

DB <- TestDb()
samples <- sample_covariate_tbl(DB) %>%
  filter(value == 'CMS4') %>%
  select(dataset, sample_id)
samples <- sample_covariate_tbl(DB) %>%
  filter(variable == 'stage' & value == 'III') %>%
  select(dataset, sample_id)
genes <- c("800", "1009", "1289", "50509", "2191", "2335", "5159")

test_that("fetch_expression limits samples correctly", {
  s.df <- collect(samples, n=Inf)

  e.sqlite <- fetch_expression(DB, samples, genes) %>% collect(n=Inf)
  e.df <- fetch_expression(DB, s.df, genes) %>% collect(n=Inf)

  ## results are same from tbl_df and tbl_sqlite `samples` parameter
  expect_equal(e.sqlite, e.df)

  ## samples limited correcly
  expect_true(setequal(paste0(e.df$dataset, e.df$sample_id),
                       paste0(s.df$dataset, s.df$sample_id)))
})

test_that("fetch_expression results converted to DGEList", {
  e <- fetch_expression(DB, samples, genes)
  y <- as.DGEList(e)
  expect_is(y, 'DGEList')

  ## check samples
  expect_is(y$samples, 'data.frame')
  expect_type(y$samples$norm.factors, 'double')
  expect_type(y$samples$lib.size, 'integer')
  expect_type(y$samples$dataset, 'character')
  expect_type(y$samples$sample_id, 'character')

  expect_is(y$genes, 'data.frame')
  expect_type(y$genes$feature_id, 'character')
  expect_type(y$genes$symbol, 'character')
})

test_that('as.DGEList assigns correct covariates', {
  e <- fetch_expression(DB, samples, genes)
  y <- as.DGEList(e, covariates=c('stage', 'sex'))

  expect_is(y, 'DGEList')
  expect_is(y$samples$sex, 'factor')
  expect_is(y$samples$stage, 'character')
})

test_that("cpm on fetch_expression result mimics cpm.DGEList", {
  e <- fetch_expression(DB, samples, genes)

  E <- cpm(e, log=FALSE, prior.count=5) %>%
    mutate(samid=paste(dataset, sample_id, sep='_')) %>%
    mcast(feature_id ~ samid, value.var='cpm')

  EL <- cpm(e, log=TRUE, prior.count=5) %>%
    mutate(samid=paste(dataset, sample_id, sep='_')) %>%
    mcast(feature_id ~ samid, value.var='cpm')

  Etbldf <- cpm(collect(e, n=Inf), log=TRUE, prior.count=5, db=DB) %>%
    mutate(samid=paste(dataset, sample_id, sep='_')) %>%
    mcast(feature_id ~ samid, value.var='cpm')

  ## counts from within database or collect should match
  expect_equal(EL, Etbldf)

  ## Test that cpm's returned with edgeR only code match
  y <- as.DGEList(e)
  YL <- cpm(y, log=TRUE, prior.count=5)
  Y <- cpm(y, log=FALSE, prior.count=5)

  E <- E[rownames(Y), colnames(Y)]
  EL <- EL[rownames(YL), colnames(YL)]

  expect_equal(E, Y)
  expect_equal(EL, YL)
})


