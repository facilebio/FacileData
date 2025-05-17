# Second take on faciledataset construction. This will enable more memory
# efficient construction of larger datasets, and more easily allow the addition
# of newer data into an already built faciledataset.
# 
# First initialize the faciledataset with the feature/assay universe that we 
# want to load it up with, then add data.
# 
# Note that this requires an initial metadata (organism)
# 
# 1. fds_initialize(directory, species):
#    a. create new directory structure
#    b. initialize sql db,
#    c. create initial hdf5 file structure
#    d. provide firsdt assay/feature map
# 2. fds_add_assay(fds, assay_info)
#    a. create assay map from feature descriptor and feature universe
# 3. fds_add_assay_data(fds, assay_name, data)
#    a. add assa daya to an assay already added in (2)
#    b. this will add sample-metadata if these are bioc data containers

#' Add an assay matrix to a pre-registered assay in an FDS.
#' 
#' This functions adds data to an already registered assay. The rownames() of
#' the assay data must already exist in the feature space that is assigned to
#' the assay.
#' 
#' Rows in `assay_matrix` that do not exist in the registered feature space
#' will be ignored.
#' 
#' @export
#' @param x a FacileDataSet
#' @param data_matrix the matrix of assay_data to add. `rownames(data_matrix)`
#'   must be already existing `feature_id`s in the registered `assay_name`
#' @param assay_name the name of a pre-registered assay in `x`
#' @param dataset_name a string indicating the `dataset` value for this data.
#'   Recall that dataset,sample_id (ie. dataset,colnames(assay_matrix)) are
#'   the primary key's of the data in `x`
#' @param ... stuff
#' @param missing_value what to fill in to the internal hdf5 assay matrix for
#'   features that are not present in `assay_matrix` but exist in the
#'   `assay_name` feature space. Default is 0
fds_add_assay_data <- function(x, assay_matrix, assay_name, dataset_name, ...,
                               missing_value = 0L) {
  checkmate::assert_class(x, "FacileDataSet")
  checkmate::assert_multi_class(assay_matrix, c("matrix", "Matrix"))
  checkmate::assert_character(rownames(assay_matrix))
  checkmate::assert_character(colnames(assay_matrix))
  checkmate::test_string(assay_name)
  checkmate::test_string(assay_type)
  
  assays_all <- assay_info(x)
  ainfo <- filter(assay_info(x), .data$assay == .env$assay_name)
  stopifnot('assay_name found' = nrow(ainfo) == 1L)
  
  if (ainfo$storage_mode == "integer") {
    missing_value <- as.integer(missing_value)
  } else {
    missing_value <- as.numeric(missing_value)
  }
  
  # feature logic
  features.fds <- x |> 
    assay_feature_info_tbl() |> 
    filter(
      .data$assay == ainfo$assay
    ) |> 
    select(feature_id, hdf5_index) |> 
    arrange(.data$hdf5_index) |> 
    collect()
  
  stopifnot(
    "hdf5_index is not sequential" = {
      all.equal(1:nrow(features.fds), features.fds$hdf5_index)
    }
  )
  
  features.data <- tibble(
    feature_id = checkmate::assert_character(rownames(assay_matrix)),
    data_index = seq(n())
  )
  
  # universe of features, ones that are registered and ones that are in the
  # assay matrix
  features.universe <- full_join(features.fds, features.data, by = "feature_id")
  features.unknown <- filter(features.universe, is.na(hdf5_index))
  features.missing <- filter(features.universe, is.na(data_index))
  
  a.go <- matrix(
    missing_value, 
    nrow = nrow(ffeatures), 
    ncol = ncol(assay_matrix),
    dimnames = list(ffeatures$feature_id, colnames(assay_matrix)))
  
  assay_matrix <- as.matrix(assay_matrix)
  if (ainfo$storage_mode == "integer" && storage.mode(assay_matrix) != "integer") {
    assay_matrix <- round(assay_matrix)
    storage.mode(assay_matrix) <- "integer"
  }
  checkmate::assert_choice()
  assay_registered <- assay_info()  
}

#' Register a new assay to a feature space
#' 
#' @export
#' @param x FacileDataSet
#' @param assay_name the name of the new assay to register
#' @param assay_type the type of assay (rnaseq, lognorm, etc.)
#' @param feature_type the identifer of one of the featurespaces in
#'   `feature_space()`
#' @examples
#' fds_register_assay(x, "counts", "rnaseq", "ensgid")
fds_register_assay <- function(
    x, 
    assay_name, 
    assay_type, 
    feature_type, ...,
    storage_mode = c("numeric", "integer")
) {
  checkmate::assert_class(x, "FacileDataSet")
  checkmate::assert_string(assay_name)
  checkmate::assert_string(assay_type)
  checkmate::assert_string(feature_type)
  storage_mode <- match.arg(storage_mode)
  
  # already exists?
  if (assay_name %in% assay_names(x)) {
    stop("Assay name already exists: ", assay_name)
  }
  
  fspace <- feature_space(x) |> 
    filter(.data$feature_type == .env$feature_type)
  
  if (nrow(fspace) != 1L) {
    stop("feature_type `", feature_type, "` not registered in dataset, ",
         "consider running fds_add_feature_space()")
  }
  
  # add to assay_info table
  ainfo <- tibble(
    assay = assay_name,
    assay_type = assay_type,
    feature_type = feature_type,
    nfeatures = fspace$n,
    storage_mode = storage_mode
  )
  
  ai <- append_facile_table(ainfo, x, "assay_info")
  
  afeatures <- fspace |> 
    transmute(
      assay = assay_name,
      feature_id,
      hdf5_index = seq(n())
    ) |> 
    append_facile_table(x, "assay_feature_info")
  
  list(
    assay_info = ainfo,
    features = afeatures
  )  
  # add map to assay_feature_info table to map feature_id <> hd5 idx
}

#' Add a new feature space to add assay data for
#' @export
#' @param x a faciledataset
#' @param feature_info a tibble of the feature universe, minimally must contain
#'   `feature_id`, `name`, `meta` columns. If `feature_type` is not a column,
#'   you can specify it from the `feature_type` parameter.
#' @param feature_type if provided, it will overvide the value in `feature_info`
#'   If not provided (default), function will fail if `feature_type` is not
#'   provided as a coloumn in `feature_info`
#' @param assays if provided, this needs to be a tibble with `assay_name`,
#'   `assay_type`, (and optional `description`) columns -- this will auto
#'   register this feature space with the enumerated assays.
#' @examples
#' afds <- an_fds()
#' fdir <- tempfile("fds__")
#' afeatures <- features(afds)
#' 
#' xfds <- fds_initialize(fdir, "human")
#' fnew <- fds_add_feature_space(xfds, afeatures)
fds_add_feature_space <- function(
    x,
    feature_info,
    description = NULL,
    feature_type = NULL,
    assays = NULL,
    ...) {
  checkmate::assert_class(x, "FacileDataSet")
  checkmate::assert_data_frame(feature_info, min.rows = 1L)
  req.cols <- c("feature_id", "name", "meta")
  checkmate::assert_subset(req.cols, colnames(feature_info))
  
  if (is.null(feature_type)) {
    feature_type <- checkmate::assert_character(feature_info[["feature_type"]])
    feature_type <- unique(feature_type)
  }
  checkmate::assert_string(feature_type)
  feature_info[["feature_type"]] <- feature_type
  
  dupes <- feature_info[duplicated(feature_info),]
  if (nrow(dupes) > 0) {
    stop("Duplicated features in feature_info")
  }
  
  if (is.null(description)) {
    description <- sprintf(
      "%d row feature-space for `%s` feature_type",
      nrow(feature_info),
      ftype
    )
  }
  checkmate::assert_string(description)
  
  if (feature_type %in% feature_types(x)) {
    stop("Feature type `", feature_type, "` already registered")
  }
  
  append_facile_feature_info(x, feature_info)
}

# Initialization ===============================================================

#' Initialize a new home for a FacileDataSet
#' @export
#' @examples
#' fdir <- tempfile("fds__")
#' xfds <- fds_initialize(fdir, "human")
fds_initialize <- function(
    directory,
    species,
    name = NULL,
    ...
) {
  dirs <- .fds_initialize_directory_structure(directory, ...)
  meta <- .fds_initialize_meta_file(dirs, species, name, ...)
  db <- .fds_initialize_db(dirs, meta, ...)
  h5 <- .fds_initialize_hdf5(dirs, meta, db, ...)
  
  FacileDataSet(
    dirs$directory, 
    data.fn = db$fn, 
    sqlite.fn = db$fn, 
    hdf5.fn = h5, 
    meta.fn = meta$fn, 
    anno.dir = dirs$annotation,
    validate_metadata = FALSE)
}

.fds_initialize_hdf5 <- function(dirs, meta, con, ...) {
  hd5.fn <- file.path(dirs$directory, "data.h5")
  rhdf5::h5createFile(hd5.fn)
  rhdf5::h5createGroup(hd5.fn, "assay")
  hd5.fn
}

.fds_initialize_db <- function(
    dirs, 
    meta, 
    ...,
    page_size = 2**12,
    cache_size = 2e5
) {
  db.fn <- file.path(dirs$directory, 'data.sqlite')
  con <- RSQLite::dbConnect(RSQLite::SQLite(), db.fn)
  RSQLite::dbExecute(con, "pragma temp_store=MEMORY;")
  RSQLite::dbExecute(con, sprintf("pragma page_size=%d", page_size))
  RSQLite::dbExecute(con, sprintf("pragma cache_size=%d;", cache_size))
  sql.fn <- system.file("extdata", "init", "faciledataset.sql",
                        package = "FacileData")
  db.sql <- sqlFromFile(sql.fn)
  executeSQL(con, db.sql)
  list(con = con, fn = db.fn)
}

.fds_initialize_meta_file <- function(
    dirs,
    species,
    name = NULL,
    meta_file = NULL,
    ...
) {
  checkmate::assert_list(dirs, names = "unique")
  if (is.null(meta_file)) {
    meta_file <- file.path(dirs$directory, "meta.yaml")
  }
  checkmate::assert_directory_exists(dirname(meta_file), "w")
  
  sinfo <- sparrow::species_info(species)
  if (is.null(name)) {
    name <- basename(dirs$directory)
  }
  checkmate::assert_string(name)
  
  meta <- list(
    name = name,
    organism = sinfo$species
  )
  
  yaml::write_yaml(meta, file = meta_file)
  attr(meta, "path") <- meta_file
  list(meta = meta, fn = meta_file)
}

.fds_initialize_directory_structure <- function(
    directory,
    annotation_dir = file.path(directory, "custom-annotation"),
    ...
) {
  if (checkmate::test_directory_exists(directory)) {
    stop("FDS directory already exists: ", directory)
  }
  pdir <- checkmate::assert_directory_exists(dirname(directory), "w")
  directory <- normalizePath(directory, mustWork = FALSE)
    
  out <- list(
    directory = directory,
    annotation = annotation_dir
  )
  for (dname in names(out)) {
    d <- out[[dname]]
    dgood <- dir.create(d)
    if (!dgood) {
      stop("Failed to create `", dname, "` facile directory: ", d)
    }
  }
  
  out$parent <- pdir
  out
}