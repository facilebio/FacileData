##' Convert +/- encoding of time to event to time + event.
##'
##' Time to event (and event labels) are encoded as a single real value in the
##' clinical data. The absoulte value of the value is the time to the event, and
##' a negative value indicates "event" and a positive value is right censored.
##'
##' @export
##' @param x the time to event
##' @return two column data.frame
parse_right_censored <- function(x, suffix=NULL) {
  x <- as.numeric(x)
  out <- data.frame(tte=abs(x), event=as.integer(x < 0))
  if (is.character(suffix)) {
    names(out) <- paste0(names(out), '_', sub('^[^a-zA-Z]', '', suffix))
  }
  out
}
