#' Summary of feature spaces in the FacileDataSet
#' 
#' @export
feature_space <- function(x, ...) {
  checkmate::assert_class(x, "FacileDataSet")
  x |> 
    feature_info_tbl() |> 
    count(feature_type) |> 
    collect()
}

#' Materializes a table with all feature information for a given assay.
#'
#' DEBUG: This logic is unnecessarily complex because I make sure to collect
#' all tables from the database as opposed to copying external tables in and
#' doing an inner_join in the database. I'm doing this becuase we are getting
#' name collisions on some of the temporary tables. We get errors like:
#'     Warning: Error in : Table pkdtpohpsu already exists.
#'
#' 
#' Retrieves feature information from the FacileDataStore for *either* a
#' particular `assay_name` or `feature_type`. By default this method returns
#' feature information for the features measured on the `default_assay()` of
#' this FacileDataStore.
#'
#' @export
#' @noRd
features.FacileDataSet <- function(x, assay_name = NULL, feature_type = NULL,
                                   feature_ids = NULL, ...) {
  null.aname <- onull.aname <- is.null(assay_name)
  null.ftype <- onull.ftype <- is.null(feature_type)
  
  # Is the user asking for feature information from the features measured on
  # a given assay, or for all features of a given feature_type.
  if (null.aname && null.ftype) {
    assay_name <- default_assay(x)
    null.aname <- FALSE
  }
  if (!xor(null.aname, null.ftype)) {
    stop("Must specify feature information for EITHER assay_name or ",
         "feature_type")
  }
  if (null.aname) {
    query_type <- "feature_type"
    query_value <- assert_choice(feature_type, feature_types(x))
  } else {
    query_type <- "assay_name"
  }
  
  if (!is.null(feature_ids)) {
    assert_character(feature_ids)
  }
  
  if (query_type == "feature_type") {
    out <- filter(feature_info_tbl(x), feature_type == query_value)
    if (!is.null(feature_ids) && length(feature_ids) > 0) {
      out <- filter(out, feature_id %in% feautre_ids)
    }
    out <- collect(out, n = Inf)
  } else {
    ftype <- assay_feature_type(x, assay_name)
    
    afinfo <- assay_feature_info_tbl(x) |>
      filter(assay == assay_name)
    
    if (!is.null(feature_ids) && length(feature_ids) > 0) {
      afinfo <- filter(afinfo, feature_id %in% feature_ids)
    }
    afinfo <- collect(afinfo, n=Inf)
    
    assay.info <- assay_info_tbl(x) |>
      select(assay, assay_type, feature_type) |>
      filter(assay == assay_name) |>
      collect(n = Inf)
    
    ## FIXME: consider materialized view for this
    out <- inner_join(afinfo, assay.info, by = "assay")
    
    ftype <- out$feature_type[1L]
    finfo <- filter(feature_info_tbl(x), feature_type == ftype)
    finfo <- collect(finfo, n = Inf)
    
    # FIXME: feature_id should be made unique to feature_type to simplify
    # e.g add GeneID: prefix for entrez
    # But, still we know out and finfo each only have one feature type now
    out <- inner_join(out, finfo, by = c("feature_type", "feature_id"))
  }
  
  if (onull.aname) {
    out <- select(out, -any_of(c("assay", "hdf5_index", "assay_type")))
  }
  
  as_facile_frame(out, x)
}
