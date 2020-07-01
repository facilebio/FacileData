context("freplace_na")

.def.categorical <- FacileData:::defaults.freplace_na$categorical

test_that("freplace_na handles factors", {
  data <- data.frame(
    a = rnorm(10),
    b = letters[1:10],
    c = factor(LETTERS[1:10]))
  data[3, 2:3] <- NA

  r1 <- freplace_na(data)
  expect_true(all(complete.cases(r1)))
  checkmate::expect_factor(
    r1$c,
    levels = c(head(LETTERS, nrow(data)), .def.categorical))
})

test_that("freplace_na errors on numerics unless given explicit replacement", {
  data <- data.frame(
    a = rnorm(10),
    b = letters[1:10],
    c = factor(LETTERS[1:10]))
  data[3, ] <- NA
  expect_error(freplace_na(data), "numerics.*number")

  r <- freplace_na(data, defaults = list(numeric = -1))
  expect_equal(r$a[3], -1)
  expect_equal(r$b[3], .def.categorical)
  expect_equal(as.character(r$c[3]), .def.categorical)
})

test_that("freplace_na handles custom values per column", {
  data <- data.frame(
    a = rnorm(10),
    b = letters[1:10],
    c = factor(LETTERS[1:10]))
  data[3, 2:3] <- NA

  r <- freplace_na(data, list(b = "bee"))
  expect_equal(r$b[3], "bee")
  expect_equal(as.character(r$c[3]), .def.categorical)
})

test_that("freplace_na ignores specified columns", {
  data <- data.frame(
    a = rnorm(10),
    b = letters[1:10],
    c = factor(LETTERS[1:10]))
  data[3, ] <- NA

  # Since a is numeric and has NA, this should error, but we explicitly ask to
  # skip the numeric column
  r <- freplace_na(data, list(b = "bee"), ignore = "a")
  checkmate::expect_scalar_na(r$a[3])
  expect_equal(r$b[3], "bee")
  expect_equal(as.character(r$c[3]), .def.categorical)
})
