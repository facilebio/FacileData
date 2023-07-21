# if (!exists("FDS")) FDS <- exampleFacileDataSet()
if (!exists("FDS")) FDS <- an_fds()

features <- some_features()
ctypes <- some_celltypes()

# A subset of these samples will not have assay data from the scRNAseq assay,
samples.sparse <- some_samples(sparse = TRUE)

# These samples have assay data from both scrnaseq and snrnaseq assay
samples.full <- some_samples(sparse = FALSE)

test_that("fetch_assay_data limits samples correctly", {
  s.df <- samples.full |> 
    collect(n = Inf) |> 
    arrange(dataset, sample_id)
  
  # fetch_assay_data over a FacileDataStore ------------------------------------
  adata.fds <- fetch_assay_data(FDS, features, samples = s.df)
  asamples.fds <- adata.fds |> 
    distinct(dataset, sample_id) |> 
    arrange(dataset, sample_id)
  expect_equal(asamples, s.df, check.attributes = FALSE)

  # fetch_assay_data from a facile_frame ---------------------------------------
  adata.ff <- fetch_assay_data(s.df, features)
  asamples.ff <- adata.ff |> 
    distinct(dataset, sample_id) |> 
    arrange(dataset, sample_id)
  expect_equal(asamples.ff, s.df, check.attributes = FALSE)
})

test_that("fetch_assay_data selects the right sample<>feature subset", {
  s.df <- samples.full |> 
    collect(n = Inf) |> 
    arrange(dataset, sample_id)
  
  expected_samples_features <- tidyr::expand_grid(
    sample_id = s.df$sample_id,
    feature_id = features$feature_id) |> 
    arrange(sample_id, feature_id)

  # fetch_assay_data over a FacileDataStore ------------------------------------
  adata.fds <- fetch_assay_data(s.df, features)
  asf <- adata.fds |> 
    select(sample_id, feature_id) |> 
    arrange(sample_id, feature_id)

  expect_equal(asf, expected_samples_features, check.attributes = FALSE)
})

test_that("fetch_assay_data can keep or drop samples without assay support", {
  # This is exercising the .fetch_assay_data function.
  # When we ask for assay data from a sample that doesn't have any data from
  # that assay, we got a long data.frame with NA values in things like
  # feature_id, ie. the first row below shouldn't happen:
  #   adata <- fetch_assay_data(samples.sparse, features)
  # dataset sample_id    assay    assay_type feature_type feature_id      feature_name value
  # <chr>   <chr>        <chr>    <chr>      <chr>        <chr>           <chr>        <int>
  # AKI     CNT.30_10034 NA       NA         NA           NA              NA              NA
  # AKI     CNT.32_10003 scrnaseq pseudobulk ensgid       ENSG00000119888 EPCAM         2736
  
  adata.drop <- expect_warning({
    samples.sparse |> 
      fetch_assay_data(features[1L,], drop_samples = TRUE, verbose = TRUE)
  }, sprintf("samples not found.*%s", features$assay[1]))
  expect_lt(nrow(adata.drop), nrow(samples.sparse))
  
  # when drop_samples = FALSE, the same samples that we sent in should be
  # the same samples that came back
  adata.keep <- expect_warning({
    samples.sparse |> 
      fetch_assay_data(features[1L,], drop_samples = FALSE, verbose = TRUE)
  }, sprintf("samples not found.*%s", features$assay[1]))
  expect_equal(nrow(adata.keep), nrow(samples.sparse))
  expect_equal(select(adata.keep, dataset, sample_id), samples.sparse,
               check.attributes = FALSE)
  
  # check that the samples that were dropped when drop_samples = TRUE are the
  # same ones with NA in the assay.
  samples_dropped <- samples(adata.drop, dropped = TRUE)
  dropped_expected <- adata.keep |> 
    filter(is.na(assay)) |> 
    select(dataset, sample_id)
  expect_equal(samples_dropped, dropped_expected, check.attributes = FALSE)
  expect_class(fds(samples_dropped), "FacileDataStore")
})

test_that("with_assay_data spreads data 'as expected'", {
  expected <- samples.full |>
    fetch_assay_data(features, normalized = TRUE) |>
    select(dataset, sample_id, feature_name, value) |>
    tidyr::spread(feature_name, value)

  result <- samples.full |>
    with_assay_data(features, normalized = TRUE) |>
    collect()
  expect_equal(result, expected, check.attributes = FALSE)
})

test_that("fetch_assay_data(..., aggregate = TRUE) provides scores", {
  scores <- samples.full |>
    fetch_assay_data(features, normalized = TRUE, aggregate = TRUE) |>
    arrange(sample_id, feature_name) |>
    select(dataset, sample_id, feature_id, symbol=feature_name, value) |>
    mutate(samid=paste(dataset, sample_id, sep="__"))

  dat <- samples.full |>
    fetch_assay_data(features, normalized = TRUE, as.matrix = TRUE)
  ewm <- sparrow::eigenWeightedMean(dat)$score[scores$samid]
  expect_equal(scores$value, unname(ewm))

  # test with_assay_data
  with.scores <- scores |>
    distinct(dataset, sample_id) |>
    with_assay_data(features, aggregate = TRUE)

  expect_equal(with.scores$aggregated, scores$value)
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
#   samples <- sample_info_tbl(tcga) |>
#     filter(dataset == 'BRCA') |>
#     collect
#
#   genes <- c(TIGIT='201633', CD274='29126')
#   rnaseq <- tcga |>
#     fetch_assay_data(genes, samples, 'rnaseq', normalized=TRUE)
#
#   ## don't have agilent data for all brca samples
#   agilent <- tcga |>
#     fetch_assay_data(genes, samples, 'agilent', normalized=TRUE)
#
# })
