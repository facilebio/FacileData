context("as.FacileDataSet")

test_that("We can get pdata metadata", {
  stopifnot(requireNamespace("Biobase", quitely = TRUE))
  stopifnot(requireNamespace("survival", quietly = TRUE))
  sinfo = data.frame(a = 1:4,
                     b = survival::Surv(1:4, c(1,1,0,1)),
                     stringsAsFactors = FALSE
  )
  rownames(sinfo) = letters[1:4]
  attr(sinfo, "label") = c(a = "a is a", b = "b is b")
  vals = matrix(1:16, ncol = 4, dimnames = list(LETTERS[1:4], letters[1:4]))
  es = Biobase::ExpressionSet(vals, Biobase::AnnotatedDataFrame(sinfo))

  expect_identical(
    FacileData::pdata_metadata(es),
    list(a = list(description = "a is a"),
         b = list(description = "b is b"))
  )
})

test_that("exampleFacileDataSet -> DGELists -> as.FacileDataSet", {
  efds <- exampleFacileDataSet()
  dsets <- sample_info_tbl(efds) %>%
    distinct(dataset) %>%
    pull(dataset)
  dlists <- sapply(dsets, function(dset) {
    y <- sample_info_tbl(efds) %>%
      filter(dataset == dset) %>%
      as.DGEList()
    y$samples <- transform(y$samples, group = NULL, samid = NULL)
    y$genes <- rename(y$genes, name = "symbol")
    colnames(y) <- sub(".*?__", "", colnames(y))
    y
  }, simplify = FALSE)

  outdir <- tempfile(pattern = "TestFacileDataSet")

  tfds <- as.FacileDataSet(dlists, outdir,
                           dataset_name = "TestFacileDataSet",
                           assay_name = "rnaseq",
                           assay_type = "rnaseq",
                           source_assay = "counts",
                           organism = organism(efds))

  # test tumor samples are equivalent
  tsamples.new <- filter_samples(tfds, sample_type == "tumor")
  tsamples.exp <- filter_samples(efds, sample_type == "tumor")
  res <- inner_join(
    mutate(tsamples.new, source = "test"),
    mutate(tsamples.exp, source = "orig"),
    by = c("dataset", "sample_id"))
  expect_equal(nrow(tsamples.new), nrow(res))
  expect_equal(nrow(tsamples.exp), nrow(res))

  # expect factor levels are the same
  stage.new <- with_sample_covariates(tsamples.new, "stage")
  stage.exp <- with_sample_covariates(tsamples.exp, "stage")
  expect_factor(stage.new[["stage"]])
  expect_equal(levels(stage.new[["stage"]]), levels(stage.exp[["stage"]]))
  stage.res <- inner_join(stage.new, stage.exp,
                          by = c("dataset", "sample_id"),
                          suffix = c(".new", ".exp"))
  expect_equal(nrow(stage.new), nrow(stage.exp))
  expect_equal(nrow(stage.new), nrow(stage.res))
  expect_equal(stage.res[["stage.new"]], stage.res[["stage.exp"]])
})
