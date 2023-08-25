#' @noRd
#' @export
has_assay.FacileDataStore <- function(x, assay_name = NULL, ...) {
  assert_facile_data_store(x)
  assert_character(assay_name, null.ok = FALSE)
  is.element(assay_name, assay_names(x))
}

#' If assay_name is NULL, we iterate over all assays to see if we got it.
#' @noRd
#' @export
has_assay.facile_frame <- function(x, assay_name= NULL,
                                   prefix = "has_", ...) {
  assert_string(prefix)
  if (is.null(assay_name)) {
    assay_name <- assay_names(fds(x))
  }
  out <- x
  for (aname in assay_name) {
    colname <- paste0(prefix, aname)
    asi <- assay_sample_info(x, aname, drop_samples = FALSE)
    asi[[colname]] <- !is.na(asi$assay) & asi$assay == aname
    res <- select(asi, dataset, sample_id, {{colname}})
    out <- left_join(out, res, by = c("dataset", "sample_id"))
  }
  out
}

#' Removes samples without specific assay support
#' 
#' Samples that do not have data from `assay_name` are dropped, and can
#' retrieved by `samples(out, dropped = TRUE)`
#'
#' @export
#' @param x a facile_frame with the same columns it was sent in with
#' @param assay_name the name of the assay
filter_by_assay_support <- function(x, assay_name, ...) {
  assert_class(x, "facile_frame")
  assert_choice(assay_name, assay_names(fds(x)))
  
  # make sure our prefix creates a unique new column name
  prefixes <- sub(assay_name, "", colnames(x), fixed = TRUE)
  prefix <- tail(make.names(c(prefixes, "..haz_"), unique = TRUE), 1L)
  
  hcol <- paste0(prefix, assay_name)
  haz <- has_assay(x, assay_name = assay_name, prefix = prefix)
  dropped <- haz |> 
    filter(!.data[[hcol]]) |> 
    distinct(dataset, sample_id)
  
  if (nrow(dropped)) {
    x <- filter(x, haz[[hcol]])
  }
  
  attr(x, "samples_dropped") <- dropped
  x
}
