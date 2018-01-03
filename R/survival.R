##' Entity-attribute-value encodings for survival data.
##'
##' @details
##' Encoding of survival data in R requires two columns, one to store
##' the time-to-event and another to indicate if there was an "event" at stored
##' time, or if it was censored. A `FacileDataSet` stores these two `pData`
##' columns into one "value" column in its entity-atribute-value
##' `sample_covariate` table.
##'
##' The `encode_right_censored` function takes the ime-to-event and censoring
##' vectors and encodes them into a single signed time-to-event numeric value.
##' Positive values indicate an event, and negative value are censored.
##'
##' The `decode_right_censored` function re-instatiaes the two-column R-native
##' storage of this data.
##'
##' @md
##' @export
##' @rdname eav-right-censor
##' @seealso [create_eav_metadata()]
##'
##' @param time `numeric` time to event
##' @param event 0/1 vector encoded in the "R sense". "1" is an event, "0" is
##'   right censored.
##' @param sas.encoding Is the 'event' vector "SAS encoded"? In the SAS world,
##'   1 means censored, and 0 is event. This is `FALSE` by default.
##' @return returns a numeric vector that combines time-to-event and censoring
##'   info (sign of the value).
encode_right_censored <- function(time, event, sas.encoding=FALSE) {
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

##' @md
##' @export
##' @rdname eav-right-censor
##' @param x the time to event
##' @param suffix adds `_<suffix>` to the `tte` and `event` columns of the
##'   outgoing `data.frame`
##' @return two column `data.frame` with `tte(_SUFFIX)?` and `event(_SUFFIX)?`
##'   columns.
decode_right_censored <- function(x, suffix=NULL, sas.encoding=FALSE) {
  x <- as.numeric(x)
  out <- data.frame(tte=abs(x))
  if (sas.encoding) {
    out[['event']] <- as.integer(x < 0)
  } else {
    out[['event']] <- as.integer(x > 0)
  }
  if (is.character(suffix)) {
    names(out) <- paste0(names(out), '_', sub('^[^a-zA-Z]', '', suffix))
  }
  out
}
