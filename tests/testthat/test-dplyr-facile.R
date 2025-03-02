test_that("group_split. returns a named version of dplyr::group_split()", {
  dsplit <- datasets::iris |> dplyr::group_split(Species)
  
  # split an ungrouped data.frame and define split in function call
  nsplit <- datasets::iris |> group_split.(Species)
  
  # split a grouped data.frame
  gsplit <- datasets::iris |> dplyr::group_by(Species) |> group_split.()
  
  expect_null(names(dsplit))
  expect_equal(names(nsplit), levels(datasets::iris$Species))
  expect_equal(names(nsplit), names(gsplit))
})

test_that("group_map. returns named version of dplyr::group_map", {
  # no names
  inrows <- datasets::iris |> group_by(Species) |> group_map(~ nrow(.x))
  # with names
  nnrows <- datasets::iris |> group_by(Species) |> group_map.(~ nrow(.x))
  
  expect_null(names(inrows))
  expect_equal(names(nnrows), levels(datasets::iris$Species))
  expect_equal(unlist(nnrows, use.names = FALSE), unlist(inrows))
})