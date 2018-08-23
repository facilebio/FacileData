library(testthat)
library(survival)

context("Coercion among Surv, cSurv and character")

test_that("We can convert among Surv, cSurv and character", {
    x = Surv(c(14,12,3), event = c(1,0,1))
    y = as(x,"cSurv")
    z = as(y, "Surv")
    x2 = as.character(x)
    z2 = as(x2, "Surv")
    expect_identical(y, structure(as.character(x), class = "cSurv"))
    expect_identical(z,x)
    expect_identical(x2, unclass(y))
    expect_identical(z2, x)
})
