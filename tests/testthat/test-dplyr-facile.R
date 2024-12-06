test_that("group_split. returns a named version of dplyr::group_split()", {
  dsplit <- iris |> dplyr::group_split(Species)
  
  # split an ungrouped data.frame and define split in function call
  nsplit <- iris |> group_split.(iris, Species)
  
  # split a grouped data.frame
  gsplit <- iris |> dplyr::group_by(Species) |> group_split.()
  
  expect_null(names(dsplit))
  expect_equal(names(nsplit), levels(iris$Species))
  expect_equal(names(nsplit), names(gsplit))
})

test_that("group_map. returns named version of dplyr::group_map", {
  # no names
  inrows <- iris |> group_by(Species) |> group_map(~ nrow(.x))
  # with names
  nnrows <- iris |> group_by(Species) |> group_map.(~ nrow(.x))
  
  expect_null(names(inrows))
  expect_equal(names(nnrows), levels(iris$Species))
  expect_equal(unlist(nnrows, use.names = FALSE), unlist(inrows))
})