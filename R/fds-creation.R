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
#' @param add_samples if TRUE (default) the dataset,sample combinations are
#'   added to the sample_info table
#' @param missing_value what to fill in to the internal hdf5 assay matrix for
#'   features that are not present in `assay_matrix` but exist in the
#'   `assay_name` feature space. Default is 0
fds_add_assay_data <- function(
    x,
    assay_matrix,
    assay_name,
    dataset_name,
    ...,
    add_samples = TRUE,
    missing_value = 0L,
    chunk_rows = 5000,
    chunk_cols = "ncol",
    chunk_compression = 4
) {
  checkmate::assert_class(x, "FacileDataSet")
  checkmate::assert_string(assay_name)
  checkmate::assert_string(dataset_name)
  checkmate::assert_flag(add_samples)
  
  # for some reason tximport creates SummarizedExperiments with `assay()` matrices
  # as actual data.frames
  checkmate::assert_multi_class(assay_matrix, c("matrix", "Matrix", "data.frame"))
  checkmate::assert_character(rownames(assay_matrix))
  checkmate::assert_character(colnames(assay_matrix))
  
  sample_id <- checkmate::assert_character(
    colnames(assay_matrix),
    len = ncol(assay_matrix),
    unique = TRUE)
  feature_id <- checkmate::assert_character(
    rownames(assay_matrix),
    len = nrow(assay_matrix),
    unique = TRUE)
  
  samples.data <- tibble(
    assay = assay_name,
    dataset = dataset_name,
    sample_id = sample_id
  )
  
  features.data <- tibble(
    feature_id = feature_id,
    data_index = seq_along(feature_id)
  )
  
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
    "hdf5_index is not sequential [this should be impossible]" = {
      all.equal(1:nrow(features.fds), features.fds$hdf5_index)
    }
  )
  
  
  # universe of features, ones that are registered and ones that are in the
  # assay matrix
  features.universe <- full_join(features.fds, features.data, by = "feature_id")
  features.unknown <- filter(features.universe, is.na(hdf5_index))
  features.missing <- filter(features.universe, is.na(data_index))
  
  if (nrow(features.unknown) > 0) {
    warning(nrow(features.unknown), " input features not registered: ignoring")
  }
  
  features.go <- features.universe |> 
    filter(!is.na(data_index), !is.na(hdf5_index))
  if (nrow(features.go) == 0) {
    stop("No features in assay_matrix match the features in registered assay")
  }
  
  assay_matrix.all <- as.matrix(assay_matrix)
  assay_matrix <- assay_matrix.all[features.go$data_index,,drop = FALSE]
  if (storage.mode(assay_matrix) == "logical") {
    storage.mode(assay_matrix) <- "integer"
  }
  checkmate::assert_choice(storage.mode(assay_matrix), c("integer", "double"))
  
  if (ainfo$storage_mode == "integer" && storage.mode(assay_matrix) != "integer") {
    assay_matrix <- round(assay_matrix)
    storage.mode(assay_matrix) <- "integer"
  }
  
  # This probably shouldn't be here, but this whole universe was built with
  # initial assumption that we really only care about NGS data
  if (ainfo$assay_type %in% c("rnaseq", "isoseq", "pseudobulk")) {
    libsize <- colSums(assay_matrix)
    normfactor <- edgeR::calcNormFactors(assay_matrix, lib.size = libsize)
  } else {
    libsize <- -1
    normfactor <- -1
  }
  
  a.go <- matrix(
    missing_value, 
    nrow = nrow(features.fds), 
    ncol = ncol(assay_matrix),
    dimnames = list(features.fds$feature_id, samples.data$sample_id)
  )
  
  if (is.character(chunk_cols) && chunk_cols == "ncol") {
    xchunk.cols <- ncol(a.go)
  } else if (is.integer(chunk_cols)) {
    xchunk.cols <- min(chunk_cols, ncol(a.go))
  } else {
    stop("Unrecognized value for chunk.cols: ", chunk_cols)
  }
  
  xchunk.rows <- min(chunk_rows = 5000, nrow(a.go))
  chunk <- c(xchunk.rows, xchunk.cols)
  dname <- sprintf('assay/%s/%s', assay_name, dataset_name)
  
  a.go[features.go$hdf5_index, ] <- assay_matrix
  
  if (chunk_compression == 0) {
    rhdf5::h5createDataset(
      hdf5fn(x),
      dname,,
      dim(a.go),
      storage.mode = ainfo$storage_mode,
      level = 0,
      filter = "NONE"
    )
  } else {
    rhdf5::h5createDataset(
      hdf5fn(x),
      dname,
      dim(a.go),
      storage.mode = ainfo$storage_mode,
      chunk = chunk,
      level = chunk_compression
    )
  }

  rhdf5::h5write(a.go, file = hdf5fn(x), name = dname)
  
  asi <- samples.data |> 
    mutate(
      hdf5_index = seq(ncol(a.go)),
      libsize = libsize,
      normfactor = normfactor
    ) |> 
    append_facile_table(x, "assay_sample_info")
  
  if (add_samples) {
    samples_added <- samples.data |> 
      distinct(dataset, sample_id) |> 
      mutate(parent_id = "") |> 
      append_facile_table(x, "sample_info")
  } else {
    samples_added <- tibble(
      dataset = character(), 
      sample_id = character(),
      parent_id = character())
  }
  
  out <- list(
    assay_sample_info = asi,
    samples_added = samples_added,
    features_unknoun = features.universe,
    reatures_missing = features.missing)
  invisible(out)
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
    feature_type, 
    description = NULL, 
    storage_mode = c("numeric", "integer"),
    ...
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
  if (!feature_type %in% feature_types(x)) {
    stop(sprintf("feature_type `%s` not registered", feature_type))
  }
  
  assay_features <- features(x, feature_type = feature_type) |> 
    arrange(feature_id) |> 
    collect() |> 
    transmute(
      assay = assay_name,
      feature_id,
      hdf5_index = seq(n())
    )
  
  if (is.null(description)) {
    description <- sprintf(
      "`%s` [%s] assay for `%s` features",
      assay_name,
      assay_type,
      feature_type
    )
  }
  checkmate::assert_string(description)
  
  # add to assay_info table
  ainfo <- tibble(
    assay = assay_name,
    assay_type = assay_type,
    feature_type = feature_type,
    description = description,
    nfeatures = nrow(assay_features),
    storage_mode = storage_mode
  )
  
  ai <- append_facile_table(ainfo, x, "assay_info")
  af <- append_facile_table(assay_features, x, "assay_feature_info")
  
  # add new subdirectory in hdf5 file
  h5name <- sprintf("assay/%s", assay_name)
  stopifnot(
    "hdf5 assay creation success" = rhdf5::h5createGroup(hdf5fn(x), h5name)
  )
  
  invisible(list(assay_info = ai, features = af))
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
    stop(nrow(dupes), " duplicated features found in feature_info")
  }
  
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