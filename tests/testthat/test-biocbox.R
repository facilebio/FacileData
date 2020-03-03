context("biocbox")

if (!exists("FDS")) FDS <- exampleFacileDataSet()

samples <- sample_covariate_tbl(FDS) %>%
  filter(variable == 'stage' & value == 'III') %>%
  select(dataset, sample_id) %>%
  collect()
genes <- local({
  out <- c("800", "1009", "1289", "50509", "2191", "2335", "5159")
  feature_info_tbl(FDS) %>%
    filter(feature_id %in% out) %>%
    collect() %>%
    pull(feature_id)
})

# boxes and their associated packages
box.info <- FacileData:::.biocboxes %>%
  select(class, package) %>%
  distinct()

test_that("fetch_assay_data results converted to biocboxes", {
  scovs <- samples %>%
    with_sample_covariates() %>%
    as.data.frame()
  rownames(scovs) <- paste(scovs$dataset, scovs$sample_id, sep = "__")

  e <- fetch_assay_data(FDS, genes, samples, as.matrix = TRUE)
  e <- e[, rownames(scovs)]

  for (i in seq(nrow(box.info))) {
    class <- box.info[["class"]][i]
    package <- box.info[["package"]][i]

    rnaseq.compat <- is.element(
      "rnaseq",
      filter(FacileData:::.biocboxes, .data$class == .env$class))

    if (rnaseq.compat) {
      bb <- biocbox(samples, features = genes, class = class)
    } else {
      bb <- expect_warning({
        biocbox(samples, features = genes, class = class)
      }, "not compatible", info = class)
    }

    expect_is(bb, class, info = class)
    checkmate::expect_set_equal(rownames(bb), genes, info = class)
    checkmate::expect_set_equal(colnames(bb), colnames(e), info = class)

    bb <- bb[rownames(e), colnames(e)]

    # Check assay data is same
    expect_equal(adata(bb), e, check.attributes = FALSE, info = class)

    # Check that sample covariates are the same
    pdat <- as.data.frame(pdata(bb))
    checkmate::expect_subset(colnames(scovs), colnames(pdat), info = class)
    expect_equal(pdat[, colnames(scovs)], scovs, info = class,
                 check.attributes = FALSE)
  }
})

test_that("biocbox appends custom covariates from input sample table", {
  custom.covs <- samples %>%
    mutate(var1 = rnorm(nrow(samples)), var2 = sample(letters, nrow(samples)))

  bb <- biocbox(custom.covs, features = genes)

  cmp <- inner_join(custom.covs, bb$samples, by = c("dataset", "sample_id"))
  expect_equal(nrow(cmp), nrow(custom.covs))
  expect_equal(cmp$var1.x, cmp$var1.y)
  expect_equal(cmp$var2.x, cmp$var2.y)
})


test_that("biocbox appends custom covariates from sample_covariates param", {
  custom.covs <- samples %>%
    mutate(var1 = rnorm(nrow(samples)), var2 = sample(letters, nrow(samples)))

  bb <- biocbox(select(custom.covs, dataset, sample_id),
                features = genes,
                sample_covariates = custom.covs)

  cmp <- inner_join(custom.covs, bb$samples, by = c("dataset", "sample_id"))
  expect_equal(nrow(cmp), nrow(custom.covs))
  expect_equal(cmp$var1.x, cmp$var1.y)
  expect_equal(cmp$var2.x, cmp$var2.y)
})
