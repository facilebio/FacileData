context("as.FacileDataSet")

es <- multiGSEA::exampleExpressionSet(do.voom=FALSE)
esl <- list(first=es, second=es)
colnames(esl[['second']]) <- paste0('two_', colnames(esl[['second']]))

test_that("yaml covariate entry entries generate successfully", {

})

test_that("single ExpressionSet converts to FacileDataSet", {

})

test_that("list of ExpressionSets convert to FacileDataSet", {

})
