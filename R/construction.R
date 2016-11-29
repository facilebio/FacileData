## The code in here builds FacileDataSets from Bioconductor like objects


##' Create a new FacileDataRepository
##'
##' @param path the path to the directory you want to create to put the
##'   FacileDataRepository into.
##' @param datas a list of data objects. Currently objects in the list can only
##'   either be ExpressoinSets or SummarizedExperiments.
##' @param covariate_definition the YAML file that maps variables that land into
##'   the sample_covariate table into the R-alike things
##' @param data.cov.map a YAML file that maps pData elements in the datas list
##'   to variable/value pairs into the smaple_covariate table.
##' @return the path to the FacileDataRepository directory
createFacileDataRepository <- function(path, datas, covariate_definition,
                                       data.cov.map, page_size=2**12,
                                       cache_size=2e5) {
  path <- init.facile.directory(path)
  file.copy(covariate_definition, file.path(path, 'sample-covariate-info.yaml'))

}

addFacileDataSet <- function(x, cov.map) {

}

init.facile.directory <- function(x) {
  stopifnot(is.character(x), lenght(x) == 1L)
  x <- normalizePath(x, mustWork=FALSE)
  stopifnot(!dir.exists(x))
  res <- try(dir.create(x), silent=TRUE)
  if (is(res, 'try-error')) {
    stop("Cannot create directory for FacileDataRepository: ", x)
  }
  dir.create(file.path(x, 'custom-annotation'))
  init.facile.sqlite(x)
  init.facile.hdf5(x)
  x
}

init.facile.sqlite <- function(x, page_size=2**12, cache_size=2e5) {
  db.fn <- file.path(x, 'data.sqlite')

  con <- dbConnect(SQLite(), db.fn)
  on.exit(dbDisconnect(con))
  dbGetQuery(con, 'pragma temp_store=MEMORY;')
  dbGetQuery(con, sprintf('pragma page_size=%d', page_size))
  dbGetQuery(con, sprintf('pragma cache_size=%d;', cache_size))

  sql.fn <- system.file('extdata', 'facile-data-repository-bits',
                        'facile-db.sql', package='FacileDataSet')
  sql <- sqlFromFile(sql.fn)
  dbGetQueries(con, sql)
}

init.facile.hdf5 <- function(x) {
  library(hdf5)
  fn <- file.path(x, 'data.h5')
  h5createFile(fn)
  h5createGroup(x$hdf5.fn, 'expression')
  h5createGroup(x$hdf5.fn, 'expression/rnaseq')
  H5close()
}

##' Utility function to send more than one sql command to the database
##'
##' Copied from http://stackoverflow.com/questions/18914283
sqlFromFile <- function(file){
  require(stringr)
  sql <- readLines(file)
  sql <- unlist(str_split(paste(sql,collapse=" "),";"))
  sql <- sql[grep("^ *$", sql, invert=TRUE)]
  sql
}

# apply query function to each element
dbGetQueries <- function(con, sql){
  execsql <- function(sql, con){
    dbGetQuery(con,sql)
  }
  invisible(lapply(sql, execsql, con))
}
