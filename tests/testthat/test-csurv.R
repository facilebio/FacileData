library(testthat)
library(survival)

context("Coercion among Surv, cSurv and character")

test_that("We can convert among Surv, cSurv and character", {
    a = Surv(c(14,12,3), event = c(1,0,1))
    b = as(a,"character")
    c = as(b, "Surv")
    expect_identical(a,c)

    d = Surv(c(14,12,3), event = c(1,0,1))
    e = as(d,"cSurv")
    f = as(e, "Surv")
    expect_identical(d,f)

    g = Surv(c(14,12,3), event = c(1,0,1))
    h = as(g,"cSurv")
    i = as(h, "character")
    expect_identical(as(i,"cSurv"), h)
    expect_identical(as(h,"Surv"), g)
    expect_identical(as(i,"Surv"), g)
    expect_identical(as(i,"cSurv"), h)

})
