# A `FacileDataSet` stores the `pData` across its samples in an
# entity-attribute-value `sample_covariate` table in the internal SQLite
# database.
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
#' @seealso [covariate_meta_info()], [covariae_defnitions()]
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
eav_decode_real <- function(x, attrname = character(), def = list(), ...) {
  out <- as.numeric(x)
  n.na <- sum(is.na(out))
  if (n.na > 0L) {
    msg <- "%d (%.2f) values in `%s` covariate failed conversion to numeric"
    warning(sprintf(msg, n.na, n.na / length(x), attrname))
  }

  out
}

eav_encode_real <- function(x, ...) {
  stopifnot(is.numeric(x))
  out <- as.character(x)
  attr(out, "eavclass") <- "real"
  out
}

eav_encode_numeric <- eav_encode_real
eav_encode_integer <- eav_encode_real

#' @rdname simple-eav-decode-functions
eav_encode_logical <- function(x, ...) {
  stopifnot(is.logical(x))
  out <- as.character(as.integer(x))
  attr(out, "eavlcass") <- "logical"
  out
}

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

#' Entity-attribute-value decoding for categorical (character) values.
#'
#' This is essentially a pass through function for categorical/character
#' values in the EAV table. If the `def` list contains a `levels` entry, then
#' the returned value is converted to a factor, with the levels in the order
#' as defined in `def$levels`. If more levels appear in `x` than exist in
#' `def$levels` they are appended to the end of the factor levels in
#' alphabetical order.
#'
#' @inheritParams eav_decode_real
#' @rdname simple-eav-decode-functions
eav_decode_categorical <- function(x, attrname=character(), def=list(), ...) {
  out <- as.character(x)
  if (is.character(def$levels)) {
    ## protect against NAing a long list of values in the event that the
    ## levels provided in the meta.yaml file don't include all levels observed
    ## here
    lvls <- unique(c(def$levels, setdiff(out, def$levels)))
    out <- factor(out, lvls)
  }
  out
}

eav_encode_categorical <- function(x, ...) {
  stopifnot(is.factor(x) || is.character(x))
  out <- as.characer(x)
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
#' columns into one "value" column in its entity-atribute-value
#' `sample_covariate` table.
#'
#' The `encode_right_censored` function takes the ime-to-event and censoring
#' vectors and encodes them into a single signed time-to-event numeric value.
#' Positive values indicate an event, and negative value are censored.
#'
#' The `decode_right_censored` function re-instatiaes the two-column R-native
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
            'on the way out',  immediate.=TRUE)
    event[isna] <- 1L
  }
  if (!all(event %in% c(0L, 1L))) {
    stop("values in 'event' that are not 0 or 1")
  }
  if (sas.encoding) {
    ## SAS is backgwards
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
#' For simple `pData` covariates, each column is treated independantly from the
#' rest. There are some types of covariates which require multiple columns for
#' proper encoding, such as encoding of survival information, which requires
#' a pair of values that indicate the "time to event" and the status of the
#' event (death or censored). In these cases, the caller needs to provide an
#' entry in the `covariate_def` list that describes which `pData` columns
#' (`varname`) goes into the single facile covariate value.
#'
#' Please refer to the **Encoding Survival Covariates** section for a more
#' detailed description of how to define encoding survival informaiton into the
#' EAV table using the `covariate_def` parameter. Further examples of how to
#' encode other complex atributes will be added as they are required, but you
#' can reference the **Encoding Arbitrarily Complex Covariates** section for
#' some more information.
#'
#' @section Encoding Survival Covariates:
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
#'     arguments=c(time="tte_OS", event="event_OS"),
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
#' Note that the `varname` entry in the `meta.yaml` file isn't used for
#' anything. It is just kept for posterity.
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
#' @param covariate_def a named list of covariate definitions. The names of
#'   this list are the names the covariates will be called in the target
#'   `FacileDataSet`. The values of the list are:
#'   * `varname`: a `character()` of the column name(s) in `x` that this
#'      sample covariate was derived from. If more than one column is to be used
#'      for the facile covariate conversion (eg. if we are encoding survival),
#'      then provide a `length() > 1` character vector with the names of the
#'      columns in `x` that were used for the encoding. If this were encoding
#'      survival this might be `c("time", "event") columns, in that order.
#'   * label: a human readable label to use for this covariate in user facing
#'      scenarios in the facileverse.
#'   * `class`: the "facile class" of the covariate. This can either be
#'     `categorical`, `real`, or `right_censored` (for survival).
#'   * `levels`: (optional) if you want a `categorical` to be treated as a
#'     factor if it isn't already encoded as such in the `pData` itself, or if
#'     you want to rearrange the factor levels.
#'   * `type`: (optional) this is used a a "grouping" level, particularly in
#'     the FacileExplorer. Not including this won't matter.
#'     TODO: talk about covariate groupings in
#'     `FacileExplorer::categoricalCovariateSelect`
#' @return a list-of-lists that encodes the `sample_covariate` section of the
#'   `meta.yaml` file for a `FacileDataSet`.
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
eav_metadata_create <- function(x, covariate_def = list()) {
  stopifnot(is.data.frame(x))
  if (is.null(covariate_def)) covariate_def <- list()
  stopifnot(is.list(covariate_def))
  if (length(covariate_def)) {
    validate_covariate_def_list(covariate_def, x)
  }

  if ("dataset" %in% colnames(x)) x[['dataset']] <- NULL
  if ("sample_id" %in% colnames(x)) x[['sample_id']] <- NULL

  # generate generic covariate definitions for all columns
  gcd <- lapply(colnames(x), function(name) eavdef_for_column(x, name))
  names(gcd) <- colnames(x)

  # remove definitions in gcd that are provided in covariate_def
  # TODO: FIXTHIS to use arguments instead of colnames
  axe <- lapply(covariate_def, '[[', 'colnames')
  axe <- unique(unlist(axe, recursive = TRUE, use.names = FALSE))
  gcd[axe] <- NULL

  out <- c(gcd, covariate_def)
  out
}

#' Generate entity-attribute-value definition for a column in a data.frame
#'
#' Creates the miniimal list-definition for a single column in a `pData`
#' `data.frame`. This function is not exported on purpose.
#'
#' @md
#'
#' @param x a `data.frame`
#' @param column the name of the column in the `x` to parse.
#' @return a generic list-of-list defintion for `x[[column]]`
eavdef_for_column <- function(x, column) {
  vals <- x[[column]]
  if (is.null(vals)) stop("Unknown column in x: ", column)

  out <- list(arguments=list(x=column), class = "categorical",
              description="no description provided")
  # encode.name <- paste0("eav_encode_", class(x)[1L])
  # encode.fn <- tryCatch(getFunction(encode.name), error = function(e) NULL)
  # if (!is.function(encode.fn)) {
  #
  # }

  # TODO: This shouldn't be an if/else thing.
  if (is.numeric(vals)) {
    out[['class']] <- "real"
  }
  if (is.logical(vals)) {
    out[['class']] <- 'logical'
  }
  if (is.factor(vals)) {
    out[['levels']] <- levels(vals)
  }
  out
}

#' Validates that a covariate defintion list reasonably describes a data.frame
#'
#' @param x a covariate definition list-of-lists
#' @param pdata a `data.frame`
validate_covariate_def_list <- function(x, pdata) {
  # this is named a list of lists
  stopifnot(is.list(x), is.character(names(x)))
  is.lists <- sapply(x, is.list)
  stopifnot(all(is.lists))
  # the names are unique
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

  # `dataset` and `sample_id` shouldn't be in here
  illegal <- intersect(c("dataset", "sample_id"), names(x))
  if (length(illegal)) {
    stop("The following variables are protected and should not be included: ",
         paste(illegal, sep = ", "))
  }

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
#' a `FacileDataSet`. This function requires the list-of-list encodings that
#' are generated from [eav_metadata_create()] to do its thing. The caller can
#' provide their own encoding via the `eav_metadata` parameter, otherwise a
#' default one will be generated.
#'
#' @md
#' @export
#'
#' @param x a wide `pData` data.frame
#' @param eav_metadata the list-of-list covariate encodings for the EAV table
#'   of the type that is generated by [eav_metadata_create()]
#' @param covariate_def passed to [eav_metadata_create()] if `eav_metadata`
#'   parameter is not provided.
#' @return a melted EAV tagble from `x`
as.EAVtable <- function(x, eav_metadata = NULL, covariate_def = list()) {
  stopifnot(is.data.frame(x))
  assert_columns(x, c("dataset", "sample_id"))

  if (is.null(eav_metadata)) {
    eav_metadata <- eav_metadata_create(x, covariate_def)
  }
  validate_covariate_def_list(eav_metadata, x)

  eav <- lapply(names(eav_metadata), function(aname) {
    eav_encode_covariate(x, eav_metadata[[vname]], aname)
  })

  dplyr::bind_rows(eav)
}

#' Encodes column(s) from `pData` into  value
#'
#' This function is not exported, and should only be called from within the
#' [as.EAVtable()] function because we rely on validity checks that are
#' hapenning there.
#'
#' @md
#'
#' @param pdata the `pData` `data.frame`
#' @param covariate_def the singe-list-definition of this covariate
#' @param vname the name of the attribute column in the eav table
#' @return a four-column `data.frame` (dataset,sample_id,variable,value)
#'   with the encoded covariate into a single `value` column.
eav_encode_covariate <- function(pdata, covariate_def, aname = "variable") {
  assert_columns(pdata, c("dataset", "sample_id"))

  # check if there is an encode function defined.
  clazz <- covariate_def[["class"]]
  if(!is.character(clazz) && length(clazz) != 1L) {
    stop("A proper `class` covariate definition is required for: ", aname)
  }

  encode.name <- paste0("eav_encode_", clazz)
  encode.fn <- getFunction(encode.name)
  if (!is.function(encode.fn)) {
    msg <- sprintf("No `%s` function defined for variable %s (%s)",
                   encode.name, aname, clazz)
    stop(msg)
  }
  args <- covariate_def[['arguments']]
  stopifnot(is.character(args), is.character(names(args)))

  ## is here a fu
  args <- as.list(args)

}

