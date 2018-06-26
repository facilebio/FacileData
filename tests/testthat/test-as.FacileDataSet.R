context("as.FacileDataSet")

#es <- multiGSEA::exampleExpressionSet(do.voom=FALSE)
#esl <- list(first=es, second=es)
#colnames(esl[['second']]) <- paste0('two_', colnames(esl[['second']]))
#
#test_that("single ExpressionSet converts to FacileDataSet", {
#})
#
#test_that("list of ExpressionSets convert to FacileDataSet", {
#})

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
    FacileData:::pdata_metadata(es),
    list(a = list(description = "a is a"),
         b = list(description = "b is b"))
  )
})
