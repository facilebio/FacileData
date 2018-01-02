##' Convert +/- encoding of time to event to time + event.
##'
##' Time to event (and event labels) are encoded as a single real value in the
##' clinical data. The absoulte value of the value is the time to the event, and
##' a positive value indicates an "event", negative is censored.
##'
##' The encoding was defined by Wei. Positive values mean that there was an
##' event at that time, negative values mean censoring (in R, that is 1 and 0,
##' respectively).
##'
##' @export
##' @param x the time to event
##' @return two column data.frame
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

##' Encodes time and event status into single numeric vector.
##'
##' A positive value for the time means that there was an event, negative
##' value is "censored".
##'
##' @export
##' @param time \code{numeric} time to event
##' @param event 0/1 vector encoded in the "R sense". "1" is an event, "0" is
##'   right censored.
##' @param sas.encoding Is the 'event' vector "SAS encoded"? In the SAS world,
##'   1 means censored, and 0 is event. This is \code{FALSE} by default.
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
