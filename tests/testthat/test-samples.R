context("Retrieving arbitrary samples")

FDS <- exampleFacileDataSet()

test_that("fetch samples returns valid sample descriptor absent assay spec", {
  s.query <- sample_info_tbl(FDS) %>%
    collect %>%
    sample_n(10)

  s.res <- fetch_samples(FDS, s.query) %>% collect
  diffs <- s.query %>%
    anti_join(s.res, by=c('dataset', 'sample_id'))

  ## these two tables should be the same
  expect_true(nrow(diffs) == 0L)
})

# test_that("fetch_samples allows filtering covariate table as if it were wide", {
#   is.tcga <- length(grep('tcga', dbfn(FDS), ignore.case=TRUE)) > 0
#   if (is.tcga) {
#     xind <- 'BLCA'
#     dats <- 'BLCA'
#     samples <- FDS %>%
#       fetch_samples(sex == 'f',
#                     indication == 'BLCA',
#                     subtype_molecular %in% c('luminal', 'basal'))
#   } else {
#     xind <- 'bladder'
#     dats <- c('pcd', 'imvigor210')
#     samples <- FDS %>%
#       fetch_samples(sex == 'f',
#                     indication == 'bladder',
#                     subtype_molecular %in% c('luminal', 'basal'))
#   }
#   asamples <- samples %>%
#     select(dataset, sample_id) %>%
#     with_sample_covariates(c('sex', 'indication', 'subtype_molecular'))
#   asamples <- droplevels(asamples)
#
#   ## If TCGA, retrieve BLCAREAD and COAD datasets
#   expect_true(setequal(asamples$dataset, dats))
#   expect_true(all(asamples$sex == 'f'))
#   expect_true(setequal(asamples$subtype_molecular, c('luminal', 'basal')))
# })

# test_that("fetch_samples throws error when subsetting with unknown covariate", {
#   asamples <- FDS %>%
#     fetch_samples(stage == 'I') %>%
#     with_sample_covariates('stage')
#   expect_true(all(asamples$stage == 'I'))
#
#   expect_error(FDS %>% fetch_samples(stages == 'I'), 'not defined: stages')
# })
