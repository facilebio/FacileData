# A `FacileDataSet` stores the `pData` across its samples in an
# entity-attribute-value `sample_covariate` table in the internal database.
#
# The functions in this file facilitate conversion of wide pdata to long EAV
# tables, and handles the value encoding and decoding into and out of the
# database and R.

# EAV -> pData utility functions ===============================================

#' Retrieve the meta information about a covariate for EAV decoding
#'
#' Mappings that define attribute-value encodings into R-native objects are
#' stored in a `FacileDataSet`'s `meta.yaml` file, in the `sample_covariate`
#' section.
#'
#' @md
#' @export
#'
#' @param covariate the name of the covariate
#' @param .fds the `FacileDataSet`
#' @param covdefs The `covariate_definitions(.fds)` list
#' @return a list of covariate information with the following elements:
#'   `$name`, `$type`, `$class`, `$description`,
#'   `$label`, `$is.factor`, (and maybe `$levels`)
covariate_meta_info <- function(covariate, .fds, covdefs=NULL) {
  if (is.null(covdefs)) {
    stopifnot(is.FacileDataSet(.fds))
    covdefs <- covariate_definitions(.fds)
  }
  assert_covariate_definitions(covdefs)
  meta <- covdefs[[covariate]]
  if (!is.list(meta)) {
    stop("Covariate `", covariate, "` not found in covariate_definition file")
  }
  meta[["name"]] <- covariate
  meta[["is.factor"]] <- meta$class == 'categorical' &&
    is.character(meta$levels) &&
    length(meta$levels) > 0L
  meta
}

#' Casts the character EAV values to their R-native defined types.
#'
#' For most things, a single value will be returned from each cast, but in the
#' case of "time_to_event" data, the value is expended to a two column
#' data.frame with a `tte_<covariate>` column for time to event, and an
#' `event_<covariate>` column to indicate event (1) or right censored (2).
#'
#' The mechanics of how values in the `sample_covariate` table are converted
#' into R objects are handled by the information stored in the
#' `FacileDataSets`'s `meta.yaml` file.
#'
#' @md
#' @export
#' @importFrom methods getFunction
#' @seealso [covariate_meta_info()], [covariate_definitions()]
#'
#' @param covariate the name of the covariate
#' @param values the covariate values (which is a `character`) as it is
#'   pulled from the database.
#' @param cov.def the un-yamled covariate definitions, if missing we rely on
#'   pulling this out from the `FacileDataSet` object `.fds`
#' @param .fds If `missing(cov.def)`, this is the `FacileDataSet` to
#'   get the covariate definitions from.
#' @return values cast to appropriate type if a valid definition was found for
#'   `covariate`, otherwise values is returned "as is". Most of the time
#'   this is a single vector, but others it can be a data.frame (for
#'   `right_censored` data, for instance)
cast_covariate <- function(covariate, values, cov.def, .fds) {
  if (missing(cov.def)) {
    stopifnot(is.FacileDataSet(.fds))
    cov.def <- covariate_definitions(.fds)
  }
  stopifnot(is(cov.def, 'CovariateDefinitions'))
  stopifnot(is.character(values))
  stopifnot(is.character(covariate) && length(covariate) == 1L)

  if (covariate != '.dummy.') {
    def <- cov.def[[covariate]]
    if (is.null(def)) {
      warning("Covariate definition not found for ", covariate,
              " casting assumed to be categorical", immediate. = TRUE)
      def <- list(class="categorical")
    }

    clazz <- def$class
    stopifnot(is.character(clazz), length(clazz) == 1L)

    decode.fn <- paste0("eav_decode_", clazz)
    dfn <- tryCatch(getFunction(decode.fn), error = function(e) NULL)
    if (!is.function(dfn)) {
      warning(decode.fn, " casting function not found, please define one.",
              " Casting currently assumed to be categorical", immediate. = TRUE)
      def$class <- 'categorical'
      dfn <- eav_decode_categorical
    }

    values <- dfn(values, covariate, def)
  }

  values
}

# Specific EAV -> pData conversion functions -----------------------------------

#' Entity-attribute-value decoding for real values.
#'
#' This is a simple function to handle converting numeric values in the EAV
#' table to numeric data in R.
#'
#' @rdname simple-eav-decode-functions
#' @param x the values column from the `EAV` table for this covariate
#' @param attrname the name of "attribute" (covariate) in the EAV table.
#' @param def the `covariate_definition` list for this covariate
#' @return a `numeric` vector of `length(x)`
#' @rdname simple-eav-decode-functions
#' @export
eav_decode_real <- function(x, attrname = character(), def = list(), ...) {
  out <- as.numeric(x)
  n.na <- sum(is.na(out))
  if (n.na > 0L) {
    msg <- "%d (%.2f) values in `%s` covariate failed conversion to numeric"
    warning(sprintf(msg, n.na, n.na / length(x), attrname))
  }

  out
}

#' @rdname simple-eav-decode-functions
#' @export
eav_encode_real <- function(x, ...) {
  stopifnot(is.numeric(x))
  out <- as.character(x)
  attr(out, "eavclass") <- "real"
  out
}

eav_encode_numeric <- eav_encode_real
eav_encode_integer <- eav_encode_real

#' @rdname simple-eav-decode-functions
#' @export
eav_encode_logical <- function(x, ...) {
  stopifnot(is.logical(x))
  out <- as.character(as.integer(x))
  attr(out, "eavlcass") <- "logical"
  out
}

#' @rdname simple-eav-decode-functions
#' @export
eav_decode_logical <- function(x, attrname = character(), def = list(), ...) {
  out <- as.logical(as.integer(x))
  isna <- is.na(out)
  nna <- sum(isna)
  if (nna > 0) {
    msg <- "%d values unsuccessfully converted to logical `%s`, kept as as NA"
    warning(sprintf(msg, nna, attrname), immediate. = TRUE)
  }
  out
}

#' @rdname simple-eav-decode-functions
#' @export
eav_encode_cSurv <- function(x, ...) {
    stopifnot(is(x, "cSurv"))
    out <- as(x, "character")
    attr(out, "eavclass") <- "cSurv"
    out
}

#' @rdname simple-eav-decode-functions
#' @export
eav_decode_cSurv <- function(x, attrname = character(), def = list(), ...) {
    as(unclass(x), "cSurv")
}

#' Entity-attribute-value decoding for categorical (character) values.
#'
#' This is essentially a pass through-function for categorical/character
#' values in the EAV table. If the `def` list contains a `levels` entry, then
#' the returned value is converted to a factor, with the levels in the order
#' as defined in `def$levels`. If more levels appear in `x` than exist in
#' `def$levels` they are appended to the end of the factor levels in
#' alphabetical order. If more levels are defined in `def$levels` than appear
#' in `x`, they are by default dropped, set `droplevels = FALSE` to keep them.
#'
#' @inheritParams eav_decode_real
#' @rdname simple-eav-decode-functions
#' @export
eav_decode_categorical <- function(x, attrname = character(), def = list(),
                                   droplevels = TRUE, ...) {
  out <- as.character(x)
  if (is.character(def[["levels"]])) {
    # protect against NAing a long list of values in the event that the
    # levels provided in the meta.yaml file don't include all levels observed
    # here
    lvls <- unique(c(def[["levels"]], setdiff(out, def[["levels"]])))
    if (droplevels) {
      lvls <- intersect(lvls, out)
    }
    out <- factor(out, lvls)
  }
  out
}

#' @rdname simple-eav-decode-functions
#' @export
eav_encode_categorical <- function(x, ...) {
  stopifnot(is.factor(x) || is.character(x))
  out <- as.character(x)
  attr(x, "eavclass") <- "categorical"
  out
}

eav_encode_character <- eav_encode_categorical
eav_encode_factor <- eav_encode_categorical

#' Entity-attribute-value encodings for survival data.
#'
#' @details
#' Encoding of survival data in R requires two columns, one to store
#' the time-to-event and another to indicate if there was an "event" at stored
#' time, or if it was censored. A `FacileDataSet` stores these two `pData`
#' columns into one "value" column in its entity-attribute-value
#' `sample_covariate` table.
#'
#' The `encode_right_censored` function takes the time-to-event and censoring
#' vectors and encodes them into a single signed time-to-event numeric value.
#' Positive values indicate an event, and negative value are censored.
#'
#' The `decode_right_censored` function re-instantiates the two-column R-native
#' storage of this data.
#'
#' @md
#' @export
#' @rdname eav-right-censor
#' @seealso [eav_metadata_create()]
#'
#' @param time `numeric` time to event
#' @param event 0/1 vector encoded in the "R sense". "1" is an event, "0" is
#'   right censored.
#' @param sas.encoding Is the 'event' vector "SAS encoded"? In the SAS world,
#'   1 means censored, and 0 is event. This is `FALSE` by default.
#' @return returns a numeric vector that combines time-to-event and censoring
#'   info (sign of the value).
eav_encode_right_censored <- function(time, event, sas.encoding=FALSE, ...) {
  event <- as.integer(event)
  isna <- is.na(event)
  if (any(isna)) {
    warning('NA values in the `event` flag, these will be set to NA again ',
            'on the way out', immediate.=TRUE)
    event[isna] <- 1L
  }
  if (!all(event %in% c(0L, 1L))) {
    stop("values in 'event' that are not 0 or 1")
  }
  if (sas.encoding) {
    ## SAS is backwards
    event <- ifelse(event == 0L, 1L, 0L)
  }
  out <- ifelse(event == 0L, -1L * time, time)
  out[isna] <- NA
  out
}

#' @md
#' @export
#' @rdname eav-right-censor
#' @param x the time to event
#' @param suffix adds `_<suffix>` to the `tte` and `event` columns of the
#'   outgoing `data.frame`
#' @param def the covariate definition for this variable
#' @return two column `data.frame` with `tte(_SUFFIX)?` and `event(_SUFFIX)?`
#'   columns.
eav_decode_right_censored <- function(x, attrname=character(), def=list(),
                                      suffix=attrname, sas.encoding=FALSE, ...) {
  x <- as.numeric(x)
  out <- data.frame(tte=abs(x))
  if (sas.encoding) {
    out[['event']] <- as.integer(x < 0)
  } else {
    out[['event']] <- as.integer(x > 0)
  }
  if (is.character(suffix) && length(suffix) == 1L) {
    names(out) <- paste0(names(out), '_', sub('^[^a-zA-Z]', '', suffix))
  }
  out
}

# pData -> EAV utility functions ===============================================

#' Create a facile covariate definition file from a sample `pData` data.frame
#'
#' Sample covariates (aka `pData`) are encoded in an
#' [entity-attribute-value (EAV) table](https://en.wikipedia.org/wiki/Entity-attribute-value_model).
#' Metadata about these covariates are stored in a `meta.yaml` file in the
#' `FacileDataSet` directory which enables the `FacileDataSet` to cast the value
#' stored in the EAV table to its native R type. This function generates the
#' list-of-list structure to represent the `sample_covariates` section of the
#' `meta.yaml` file.
#'
#' For simple `pData` covariates, each column is treated independently from the
#' rest. There are some types of covariates which require multiple columns for
#' proper encoding, such as encoding of survival information, which requires
#' a pair of values that indicate the "time to event" and the status of the
#' event (death or censored). In these cases, the caller needs to provide an
#' entry in the `covariate_def` list that describes which `pData` columns
#' (`varname`) goes into the single facile covariate value.
#'
#' Please refer to the **Encoding Survival Covariates** section for a more
#' detailed description of how to define encoding survival information into the
#' EAV table using the `covariate_def` parameter. Further examples of how to
#' encode other complex attributes will be added as they are required, but you
#' can reference the **Encoding Arbitrarily Complex Covariates** section for
#' some more information.
#'
#' @section Encoding Survival Covariates:
#'
#' UPDATE: FacileData can now use survival data encoded as a `survival::Surv` object
#' stored as a pData column. Read on for the original encoding strategy, which
#' is still implemented.
#'
#' Survival data in R is typically encoded by two vectors. One vector that
#' indicates the "time to event" (tte), and a second to indicate whether or not
#' the denoted tte is an "event" (1) or "censored" (0).
#'
#' Normally these vectors appear as two columns in an experiment's `pData`,
#' and therefore need to be encoded into the `FacileDataSet`'s EAV table. To do
#' so, the pair of vectors are turned into a signed numeric value. The absolute
#' value of the numeric indicates the "time to event" and the sign of the value
#' indicates its censoring status.
#'
#' Let's assume we have `tte_OS` and `event_OS` column that are used to encode
#' a patient's overall survival (time and censor status). To store this as an
#' "OS" covariate in the EAV table, a `covariate_def` list-of-list definition
#' that captures this encoding would look like this:
#'
#' ```
#' covariate_def <- list(
#'   OS=list(
#'     class="right_censored",
#'     arguments=list(time="tte_OS", event="event_OS"),
#'     label="Overall Survival",
#'     type="clinical",
#'     description="Overall survival in days"))
#' ```
#'
#' Note how the name of the list-entry in `covariate_def` defines the name of
#' the covariate in the `FacileDataSet`. The `class` entry for the `OS`
#' definition indicates the type of variable this is. The `arguments` section
#' is only used when encoding a wide `pData` into the EAV value column.
#' `names(arguments)` correspond to the parameters in the
#' `[eav_encode_right_censored()]` function, and their values are the columns in
#' the target `pData` that populate the respective parameters in the function
#' call. The analagous `meta.yaml` entry in the `sample_covariates` section for
#' the `"OS"` `covariate_def` entry looks like so:
#'
#' ```
#' sample_covariates:
#'   OS:
#'     class: right_censored
#'     arguments:
#'       time: tte_OS
#'       event: event_OS
#'     label: "Overall Survival"
#'     type: "clinical"
#'     description: "Overall survival in days"
#' ```
#'
#' @section Encoding Arbitrarily Complex Covariates:
#'
#' To encode a new type of complex covariate from a wide `pData` data.frame,
#' we need to:
#'
#' 1. Specify a new `class` (like `"right_censored"`) for use within a
#'    `FacileDataSet`.
#' 2. Define an `eav_encode_<class>(arg1, arg2, ...)` function which takes the
#'    R data vectors (arg1, arg2) and converts them into a single value for the
#'    EAV table.
#' 3. Define a `eav_decode_<class>(x, attrname, def, ...)` function which takes
#'    the single value in the EAV table and casts it back into the R-naive data
#'    vector(s).
#'     - `x` is the vector of (character) values from the EAV table
#'     - `attrname` is the name of the covariate in the EAV table
#'     - `def` is the definition-list for this covariate.
#'     - `...` allows each decode function to be further customized.
#'
#' @md
#' @rdname eav-metadata
#' @export
#'
#' @param x a `pData` `data.frame`
#' @param ignore the columns in `x` to not create covariate definitions for.
#'   This defaults to `c("dataset", "sample_id")` since we are in the
#'   facileverse.
#' @param covariate_def a named list of covariate definitions. The names of
#'   this list are the names the covariates will be called in the target
#'   `FacileDataSet`. The values of the list are:
#'   * `varname`: a `character()` of the column name(s) in `x` that this
#'      sample covariate was derived from. If more than one column is to be used
#'      for the facile covariate conversion (e.g. if we are encoding survival),
#'      then provide a `length() > 1` character vector with the names of the
#'      columns in `x` that were used for the encoding. If this were encoding
#'      survival this might be `c("time", "event")` columns, in that order.
#'   * `label`: a human readable label to use for this covariate in user facing
#'      scenarios in the facileverse.
#'   * `class`: the "facile class" of the covariate. This can either be
#'     `categorical`, `real`, or `right_censored` (for survival).
#'   * `levels`: (optional) if you want a `categorical` to be treated as a
#'     factor if it isn't already encoded as such in the `pData` itself, or if
#'     you want to rearrange the factor levels.
#'   * `type`: (optional) this is used a a "grouping" level, particularly in
#'     the FacileExplorer.
#' @return a list-of-lists that encodes the `sample_covariate` section of the
#'   `meta.yaml` file for a `FacileDataSet`. Each list element will have the
#'   following elements:
#'
#'   1. arguments: the name(s) of the columns from `x` used in this covariate
#'      description.
#'   2. class: `"real"`, `"categorical"`, (survival needs a bity of work)
#'   3. description: a string with minimal description
#'   4. type: this isn't really used in the dataset, but another application
#'      might want to group covariates by type.
#'
#' @examples
#' # covariate_def definition to take tte_OS and tte_event columns and turn
#' # into a facile "OS" right_censored survival covariate
#' cc <- list(
#'   OS=list(
#'     arguments=list(time="tte_OS", event="event_OS"),
#'     label="Overall Survival",
#'     class="right_censored",
#'     type="clinical",
#'     description="Overall survival in days"))
eav_metadata_create <- function(x, ignore = c("dataset", "sample_id"),
                                covariate_def = list()) {
  stopifnot(is.data.frame(x))
  if (is.null(covariate_def)) covariate_def <- list()
  stopifnot(is.list(covariate_def))
  if (length(covariate_def)) {
    # validate_covariate_def_list(covariate_def, x)
    assert_covariate_definitions(covariate_def)
  }

  cnames <- setdiff(colnames(x), ignore)

  # generate generic covariate definitions for all columns in x
  gcd <- sapply(cnames, function(name) {
    eavdef_for_column(x[[name]], name)
  }, simplify = FALSE)

  eav_metadata_merge(gcd, covariate_def)
}

#' Merge inferred and explicit covariate column metadata.
#'
#' Takes a list of (perhaps) default sets of entity-attribute metadata, as would
#' be generated from `eav_metadata_create(pData(eSet), covariate_def = NULL)`,
#' and pulls out the sister custom-definitions from the `covariate_def`
#' attribute definition list.
#'
#' If the `covariate_def` list-of-lists has information for variables not
#' found in `default_def`, ie. the definitions returned from
#' `setdiff(names(covariate_def), names(default_def))` will be added to
#' the object returned from this funciton.
#'
#' @export
#' @param default_def A list of covariate-definition-lists, as would be returned
#'   from [eav_metadata_create()]
#' @param covariate_def list of additional covariate info, such as 'label'. The
#'   variables defined here (defined by `names(covarate_def)`) need not be
#'   identical to `names(defeault_def)`.
#' @return list of column metadata lists
eav_metadata_merge <- function(default_def, covariate_def = list()) {
  # assert_covariate_definitions(default_def)

  if (is.null(covariate_def) || length(covariate_def) == 0L) {
    return(default_def)
  }
  assert_covariate_definitions(covariate_def)

  # use any covariate definitions found in `covariate_def` to override the
  # default values provided in default_def
  for (cname in names(covariate_def)) {
    def <- c(covariate_def[[cname]], default_def[[cname]])
    default_def[[cname]] <- def[!duplicated(names(def))]
  }

  # Run through the custom covariate_def list to identify if any covariates
  # are multi-column compounded elements (like survival). If so, then those
  # top-level covariate definitions are removed from the outgoing definitions.
  multi.cols <- lapply(default_def, function(def) {
    args <- def[["arguments"]]
    if (length(args) > 1) args else NULL
  })
  multi.cols <- unlist(multi.cols, use.names = FALSE)
  if (length(multi.cols)) {
    multi.dup <- multi.cols[duplicated(multi.cols)]
    if (length(multi.dup)) {
      stop("There are single columns that are used in > 1 multi-column ",
           "covariate definintions: ",
           paste(multi.dup, collapse = ","))
    }
  }
  for (mres in multi.cols) {
    if (!is.null(mres)) default_def[mres] <- NULL
  }

  default_def
}

#' Generate entity-attribute-value definition for a column in a data.frame
#'
#' Creates the minimal list-definition for a single column in a `pData`
#' `data.frame`. This function is not exported on purpose. Column descriptions
#' will be taken from the "label" attribute of data.frames or the "metadata" list
#' for DataFrames.
#'
#' @md
#'
#' @param column a vector, e.g. a column out of a pdata
#' @param column_name single character, name of the colum
#' @return a generic list-of-list definition column
eavdef_for_column <- function(column, column_name) {
  out <- list(
      arguments = list(x = column_name),
      class = "categorical",
      label = column_name,
      description = column_name,
      type = "general"
  )

  if (is.numeric(column)) {
    out[['class']] <- "real"
  }
  if (is.logical(column)) {
    out[['class']] <- 'logical'
  }
  if (is.factor(column)) {
    out[['levels']] <- levels(column)
  }
  if (is(column, "cSurv")) {
      out[['class']] <- 'cSurv'
  }
  out
}

#' Validates that a covariate defintion list reasonably describes a data.frame.
#'
#' The covariates defined in `x` must be a subset of the columns in `pdata`.
#' This method will throw an error if there is a covariate in `x` that does
#' not have a matching column in `pdata`.
#'
#' This function does not check if all columns in `pdata` have definitions in
#' `x`.
#'
#' @param x a covariate definition list-of-lists
#' @param pdata a `data.frame`
validate_covariate_def_list <- function(x, pdata) {
  # stop("Not sure about this function")
  # this is named a list of lists

  stopifnot(is.list(x), is.character(names(x)))
  is.lists <- sapply(x, is.list)
  stopifnot(all(is.lists))
  # the names of `x` are unique
  stopifnot(length(unique(names(x))) == length(x))

  # each list item has the following elements
  required.elements <- c("arguments", "class")
  # 'label', 'type' (group), and 'description' are also nice to have, but not
  # strictly required.
  has.elems <- sapply(x, function(el) required.elements %in% names(el))
  rownames(has.elems) <- required.elements
  for (wut in required.elements) {
    missed <- which(!has.elems[wut,])
    if (length(missed)) {
      msg <- sprintf("The following columns are missing '%s'\n:    %s",
                     wut, paste(colnames(has.elems)[missed], collapse=","))
      stop(msg)
    }
  }

  # # `dataset` and `sample_id` shouldn't be in here
  # illegal <- intersect(c("dataset", "sample_id"), names(x))
  # if (length(illegal)) {
  #   stop("The following variables are protected and should not be included: ",
  #        paste(illegal, sep = ", "))
  # }

  # varname entries are valid columns in `pdata`
  for (element in names(x)) {
    args <- x[[element]]$arguments
    cnames <- unlist(args)
    if (!is.character(cnames)) {
      stop(sprintf("%s:varname is not a character", element))
    }
    assert_columns(pdata, cnames)
  }

  TRUE
}

#' Convert a `pData` data.frame to a melted EAV table
#'
#' Transforms a wide `pData` data.frame into a melted EAV table for use in
#' a `FacileDataSet`. This function will also produce the list-of-list encodings
#' that are generated from [eav_metadata_create()] to do its thing as an
#' attribute of the returned object.
#'
#' If you want to provide custom definitions for the covariates in the EAVtable
#' that are different than the ones generated in [eav_metadata_create()], then
#' provie that definition list in the `covariate_def` parameter.
#'
#' @md
#' @export
#'
#' @param x a wide `pData` data.frame
#' @param covariate_def passed to [eav_metadata_create()] that is used to
#'   override default covariate definitions extracted from the columns of `x`
#' @return a melted EAV table from `x`
as.EAVtable <- function(x, ignore = c("dataset", "sample_id"),
                        covariate_def = list(), na.rm = TRUE) {
  stopifnot(is.data.frame(x))
  if (is.null(ignore)) {
    ignore <- character()
  }
  assert_subset(ignore, colnames(x))

  meta <- as.data.frame(select(x, ignore), stringsAsFactors = FALSE)
  dat <- as.data.frame(select(x, -!!ignore), stringsAsFactors = FALSE)

  # special casing "Surv" object play, for now
  for (cname in colnames(dat)) {
    vals <- dat[[cname]]
    if (is(vals, "Surv")) dat[[cname]] <- as_cSurv(vals)
  }

  eav_metadata <- eav_metadata_create(dat, ignore = NULL,
                                      covariate_def = covariate_def)
  validate_covariate_def_list(eav_metadata, dat)

  encoded <- sapply(names(eav_metadata), function(cname) {
    type <- eav_metadata[[cname]][["type"]]
    if (is.null(type)) type <- "unspecified"
    e <- eav_encode(dat[[cname]], eav_metadata[[cname]], cname)
    e[["type"]] <- type
    bind_cols(meta, e)
  }, simplify = FALSE)

  out <- as.tbl(bind_rows(encoded))
  if (na.rm) {
    out <- filter(out, !is.na(value))
  }
  attr(out, "covariate_def") <- eav_metadata
  out
}

#' Encodes column(s) from `pData` into character values
#'
#' This function is not exported, and should only be called from within the
#' [as.EAVtable()] function because we rely on validity checks that are
#' happening there.
#'
#' @md
#'
#' @param dat the vector to values to encode into an EAV table
#' @param covariate_def the single-list-definition of this covariate
#' @param vname the name of the attribute column in the eav table
#' @return a four-column `data.frame` (dataset,sample_id,variable,value)
#'   with the encoded covariate into a single `value` column.
eav_encode <- function(dat, covariate_def, varname) {
  assert_list(covariate_def)
  assert_subset(c("class"), names(covariate_def))
  assert_string(varname, min.chars = 1)

  # check if there is an encode function defined.
  clazz <- covariate_def[["class"]]
  if (!test_string(clazz)) {
    stop("A proper `class` covariate definition is required for: ", varname)
  }

  encode.name <- paste0("eav_encode_", clazz)
  encode.fn <- getFunction(encode.name)
  if (!is.function(encode.fn)) {
    msg <- sprintf(
      "No `%s` function defined for variable %s (%s)",
      encode.name, varname, clazz
    )
    stop(msg)
  }

  ## FIXME: we need to do something with args for right_censored clazz
  #args <- covariate_def[['arguments']] # Let us assume is a list of character(1)s for now
  tibble(variable = varname,
         value = encode.fn(dat),
         class = clazz)
}

#' Encodes column(s) from `pData` into character values
#'
#' This function is not exported, and should only be called from within the
#' [as.EAVtable()] function because we rely on validity checks that are
#' happening there.
#'
#' @md
#'
#' @param pdata the `pData` `data.frame`
#' @param covariate_def the single-list-definition of this covariate
#' @param vname the name of the attribute column in the eav table
#' @return a four-column `data.frame` (dataset,sample_id,variable,value)
#'   with the encoded covariate into a single `value` column.
# eav_encode_covariate <- function(pdata, covariate_def, aname = "variable") {
eav_encode_covariate <- function(dat, covariate_def, aname = "variable") {
  assert_columns(pdata, c("dataset", "sample_id"))

  # check if there is an encode function defined.
  clazz <- covariate_def[["class"]]
  if(!is.character(clazz) && length(clazz) != 1L) {
    stop("A proper `class` covariate definition is required for: ", aname)
  }

  encode.name <- paste0("eav_encode_", clazz)
  encode.fn <- getFunction(encode.name)
  if (!is.function(encode.fn)) {
      msg <- sprintf(
          "No `%s` function defined for variable %s (%s)",
          encode.name, aname, clazz
      )
    stop(msg)
  }
  ## FIXME: we need to do something with args for right_censored clazz
  #args <- covariate_def[['arguments']] # Let us assume is a list of character(1)s for now
  encode.fn(pdata[[aname]])
}
