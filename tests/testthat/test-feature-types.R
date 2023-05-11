context("Features Types")

test_that("Different classes of identifiers guessed correctly", {
  expected <- tribble(
    ~id_type,      ~id,                             ~organism,
    "refseq",      "NC_000023.11",                  "unknown",
    "refseq",      "NC_000023.10",                  "unknown",
    "refseq",      "NM_001306206.1",                "unknown",
    "refseq",      "NP_001293135.1",                "unknown",
    "refseq",      "NC_000023",                     "unknown",
    "refseq",      "NM_001306206",                  "unknown",
    "ensgid",      "ENSG00000101811",               "Homo sapiens",
    "ensgid",      "ENSMUSG00000030088",            "Mus musculus",
    "ensgid",      "ENSMUSG00000030088.2",          "Mus musculus",
    "enstid",      "ENST00000415585.6",             "Homo sapiens",
    "enstid",      "ENSMUST00000113287.7",          "Mus musculus",
    "enstid",      "ENSMUST00000113287",            "Mus musculus",
    "entrez",      "85007",                         "unknown")

  res <- infer_feature_type(expected$id)
  expect_equal(res$id, expected$id)
  expect_equal(res$id_type, expected$id_type)
  # expect_equal(res$source_organism, expected$organism)
})
