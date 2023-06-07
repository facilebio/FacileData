# if (!exists("FDS")) FDS <- exampleFacileDataSet()
if (!exists("FDS")) FDS <- testFacileDataSet()

# These were identified because they score highly on an ANOVA across celltype
features <- tibble::tribble(
  ~feature_id,        ~name,
  "ENSG00000119888",  "EPCAM",
  "ENSG00000149564",  "ESAM",
  "ENSG00000147676",  "MAL2",
  "ENSG00000128567",  "PODXL") |> 
  mutate(assay = "scrnaseq")

# These genes have strong/variable expression across these celltypes
ctypes <- c(
  "B", "C-PC", "CNT", "DCT", "PT", "TAL",
  "EC",
  "MPC",
  "POD")

# A subset of these samples will not have assay data from the scRNAseq assay,
# and we need to handle this gracefully.
samples.sparse <- FDS |> 
  filter_samples(cell_abbrev %in% ctypes)

# Manual inspection shows that the selection criteria below provides a list
# of samples that all have data from the scrnaseq assay
samples.full <- FDS |> 
  filter_samples(cell_abbrev %in% c("PT", "IMM")) |> 
  with_sample_covariates(c("cell_abbrev", "donor_id")) |> 
  filter(!donor_id %in% c("29-10008"))


test_that("fetch_assay_data limits samples correctly", {
  s.df <- collect(samples, n = Inf)
  
  e.sqlite <- fetch_assay_data(FDS, features, samples) |> collect(n=Inf)
  e.df <- fetch_assay_data(FDS, features, s.df) |> collect(n=Inf)

  ## results are same from tbl_df and tbl_sqlite `samples` parameter
  expect_equal(e.sqlite, e.df)

  ## samples limited correcly
  expect_true(setequal(paste0(e.df$dataset, e.df$sample_id),
                       paste0(s.df$dataset, s.df$sample_id)))

})

test_that("sparse sample<>assay: retention of samples without assay support", {
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
    fetch_assay_data(samples.sparse, features, drop_samples = TRUE, verbose = TRUE)
  }, sprintf("samples not found.*%s", features$assay[1]))
  expect_equal(sum(is.na(adata.drop$feature_id)), 0)
  
  adata.keep <- expect_warning({
    fetch_assay_data(samples.sparse, features, drop_samples = FALSE, verbose = TRUE)
  }, sprintf("samples not found.*%s", features$assay[1]))
  expect_gt(sum(is.na(adata.keep$feature_id)), 0)
  
})


test_that("spreading data works with_assay_data", {
  expected <- FDS |>
    fetch_assay_data(features, samples, normalized = TRUE) |>
    select(dataset, sample_id, feature_name, value) |>
    tidyr::spread(feature_name, value)
  result <- samples |>
    with_assay_data(features, normalized = TRUE, .fds = FDS) |>
    collect()
  expect_equal(result, expected, check.attributes = FALSE)
})

test_that("fetch_assay_data(..., aggregate = TRUE) provides scores", {
  scores <- FDS |>
    fetch_assay_data(features, samples, normalized = TRUE, aggregate = TRUE) |>
    arrange(sample_id, feature_name) |>
    select(dataset, sample_id, feature_id, symbol=feature_name, value) |>
    mutate(samid=paste(dataset, sample_id, sep="__"))

  dat <- FDS |>
    fetch_assay_data(features, samples, normalized = TRUE, as.matrix = TRUE)
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
