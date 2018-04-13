#' Utility function to send more than one sql command to the database
#'
#' Copied from http://stackoverflow.com/questions/18914283
#'
#' @export
#' @param file single character, name of file with SQL statements
sqlFromFile <- function(file){
  requireNamespace("stringr") || stop("Failed to require stringr")
  sql <- readLines(file)
  sql <- gsub("--.*$", '', sql) ## remove comments
  sql <- unlist(str_split(paste(sql,collapse=" "),";"))
  sql <- sql[grep("^ *$", sql, invert=TRUE)]
  sql
}

#' Execute multiple queries against a database
#' @export
#' @param con database handle
#' @param sql list of charvecs (SQL statements)
dbGetQueries <- function(con, sql){
  execsql <- function(sql, con) {
    # message(sql)
    dbGetQuery(con,sql)
  }
  invisible(lapply(sql, execsql, con))
}
