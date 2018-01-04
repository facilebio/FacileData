#' Utility function to send more than one sql command to the database
#'
#' Copied from http://stackoverflow.com/questions/18914283
#'
#' @export
sqlFromFile <- function(file){
  require(stringr)
  sql <- readLines(file)
  sql <- gsub("--.*$", '', sql) ## remove comments
  sql <- unlist(str_split(paste(sql,collapse=" "),";"))
  sql <- sql[grep("^ *$", sql, invert=TRUE)]
  sql
}

#' Execute multiple queries against a database
#' @export
dbGetQueries <- function(con, sql){
  execsql <- function(sql, con) {
    # message(sql)
    dbGetQuery(con,sql)
  }
  invisible(lapply(sql, execsql, con))
}
