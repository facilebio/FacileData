##' Creates the database
##'
##' @export
##' @importFrom edgeR cpm calcNormFactors
##' @importClassesFrom edgeR DGEList
##' @importFrom DBI dbDriver
##' @importFrom RSQLite dbConnect dbWriteTable dbSendQuery
createWarehouse <- function(db.path, datasets, sample.meta,
                            pragma.page_size=2**12, ...) {
  if (FALSE) {
    ## Run build-db.R in
    ## ~/workspace/projects/CIT/clinx/clinXwarehouse/inst/build-db/2016-04-20
    ## to load up db.fn and spiffy
    db.path <- db.fn
    datasets <- spiffy
    sample.meta <- smeta
  }

  ## Check arguments -----------------------------------------------------------
  assertPathForOutput(db.path)
  if (file.exists(db.path)) {
    stop("SQLite DB file already exists")
  }

  ## Check integrety of `datasets`
  stopifnot(is.list(datasets))
  stopifnot(all(sapply(datasets, is, 'ExpressionSet')))
  gi <- fData(datasets[[1L]])
  stopifnot(all(sapply(datasets, function(x) all.equal(fData(x), gi))))
  rownames(gi) <- NULL

  ## Setup database skeleton ---------------------------------------------------
  ## Configure DB connection
  sqlite <- DBI::dbDriver('SQLite')
  db <- RSQLite::dbConnect(sqlite, db.path)
  RSQLite::dbGetQuery(db, sprintf('pragma page_size=%d', pragma.page_size))
  RSQLite::dbSendQuery(out$con, 'pragma temp_store=MEMORY;')
  RSQLite::dbSendQuery(out$con, 'pragme cache_size=20000;')

  ## RSQLite::dbGetQuery(db, 'pragma page_size=4096')
  createGeneTable(db)
  createExpressionTables(db)
  createSampleCovariateTable(db)

  ## Populate tables -----------------------------------------------------------
  ## 1. Gene Info
  RSQLite::dbWriteTable(db, 'gene_info', gi, append=TRUE)
  RSQLite::dbSendQuery(db, 'CREATE INDEX gene_symbol ON gene_info (symbol);')
  RSQLite::dbSendQuery(db, 'CREATE INDEX gene_feature_type ON gene_info (feature_type);')

  ## 2. Expression
  initializeWithExpressionData(db, datasets)
  ## 3. sample meta (grouping)
  ## Add "now" unix time stame
  sample.meta <- sample.meta %>%
    transform(date_entered=as.numeric(Sys.time())) %>%
    as.data.frame
  RSQLite::dbWriteTable(db, 'sample_covariate', sample.meta, append=TRUE)
  RSQLite::dbSendQuery(db, 'CREATE INDEX sample_cov_var ON sample_covariate (variable);')
  RSQLite::dbSendQuery(db, 'CREATE INDEX sample_cov_val ON sample_covariate (value);')

  dbDisconnect(db)
  db.path
}

## Utility functions to update warehouse with new data -------------------------

## Things we might want
##   * addExpressionSet
##     This should add the expression data, then update the sample_stats table
##     with new normfactor estimates
##   * updateSampleMetaTable
##     just add new rows, and when we return metadata, only return the latest
##     version of the data? -- maybe

## Create and initialize the warehouse database --------------------------------

##' Creates the expression and sample_stats table, initialize with all data.
##' @importFrom reshape2 melt
##' @param db the db connection
##' @param datasets a list of expression sets to load up databases with
initializeWithExpressionData <- function(db, datasets) {
  counts <- lapply(datasets, exprs)
  counts <- do.call(cbind, counts)
  Y <- edgeR::calcNormFactors(DGEList(counts, genes=fData(datasets[[1]])))
  ## thresh <- cpm(10, mean(Y$samples$lib.size))[1L]
  ## keep <- rowSums(cpm(Y) >= thresh) >= 0.025 * ncol(Y)
  ## sum(Y$counts)
  dataset <- lapply(names(datasets), function(name) {
    rep(name, ncol(datasets[[name]]))
  })

  sample.stats <- Y$samples %>%
    mutate(dataset=unlist(unname(dataset)), sample_id=colnames(Y)) %>%
    rename(libsize=lib.size, normfactor=norm.factors) %>%
    select(dataset, sample_id, libsize, normfactor) %>%
    as.data.frame
  RSQLite::dbWriteTable(db, 'sample_stats', sample.stats, append=TRUE)

  count.df <- lapply(names(datasets), function(name) {
    reshape2::melt(exprs(datasets[[name]])) %>%
      setNames(c('feature_id', 'sample_id', 'count')) %>%
      mutate(dataset=name) %>%
      select(dataset, sample_id, feature_id, count) %>%
      mutate(sample_id=as.character(sample_id), feature_id=as.character(feature_id))
  }) %>% bind_rows %>% as.data.frame

  ## write out data.frame to `expression` table
  ## takes about 60 seconds for ~ 775 samples
  RSQLite::dbWriteTable(db, 'expression', count.df, append=TRUE)

  ## Add index on feature_id (takes about a minue, too)
  RSQLite::dbSendQuery(db, 'CREATE INDEX expression_gene ON expression (feature_id);')
  ## Add UNIQUE INDEX on what should be the PRIMARY KEY
  RSQLite::dbSendQuery(db, 'CREATE UNIQUE INDEX expression_sample_gene ON expression (dataset, sample_id, feature_id);')
  invisible(NULL)
}

createGeneTable <- function(db) {
  table.sql <- paste0(
    "CREATE TABLE gene_info (",
    "feature_id TEXT, feature_type TEXT, symbol TEXT, n_exons INTEGER, length INTEGER, source TEXT,",
    "PRIMARY KEY (feature_id));")
  RSQLite::dbSendQuery(db, table.sql)
}

##' Creates the table that holds gene counts and meta information about data
##'
##' @importFrom reshape2 melt
createExpressionTables <- function(db) {
  ## table.sql <- paste0(
  ##   "CREATE TABLE expression (",
  ##   "dataset TEXT, sample_id TEXT, feature_id TEXT, count INTEGER, ",
  ##   "PRIMARY KEY (dataset, sample_id, feature_id))")
  table.sql <- paste0(
    "CREATE TABLE expression (",
    "dataset TEXT, sample_id TEXT, feature_id TEXT, count INTEGER)")
  dbSendQuery(db, table.sql)

  table.sql <- paste0(
    "CREATE TABLE sample_stats(",
    "dataset TEXT, sample_id TEXT, libsize INTEGER, normfactor REAL, ",
    "PRIMARY KEY (dataset, sample_id))")
  dbSendQuery(db, table.sql)
}

createSampleCovariateTable <- function(db) {
  table.sql <- paste0(
    "CREATE TABLE sample_covariate(",
    "dataset TEXT, sample_id TEXT, ",
    "variable TEXT, value TEXT, ",
    "class TEXT, ", ## class: tumor_classification, clinical, response
    "type TEXT, ",  ## type real/categorical,
    "date_entered INTEGER, ", ## date stored as seconds after unix epoch
    "PRIMARY KEY (dataset, sample_id, variable))"
  )
  dbSendQuery(db, table.sql)
}
