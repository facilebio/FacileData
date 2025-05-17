#' Utility function to send more than one sql command to the database
#'
#' Copied from http://stackoverflow.com/questions/18914283
#'
#' @param file single character, name of file with SQL statements
sqlFromFile <- function(file){
  # requireNamespace("stringr") || stop("Failed to require stringr")
  sql <- readLines(file)
  sql <- gsub("--.*$", "", sql) ## remove comments
  sql <- unlist(strsplit(paste(sql, collapse = " "), ";"))
  sql <- sql[grep("^ *$", sql, invert=TRUE)]
  sql
}

#' Execute multiple queries against a database
#'
#' @param con database handle
#' @param sql list of charvecs (SQL statements)
executeSQL <- function(con, sql){
  execsql <- function(sql, con) {
    DBI::dbExecute(con, sql)
  }
  invisible(lapply(sql, execsql, con))
}
