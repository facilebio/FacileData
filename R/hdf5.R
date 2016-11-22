##' Creates an HDF5 enabled FacileDataSet from a completely sqlite-based one.
##'
##' This is a convenience function that will be axed eventually.
createExpressionHDF5 <- function(x, drop.expr.table=TRUE) {
  library(RSQLite)
  library(rhdf5)
  stopifnot(is.FacileDataSet(x))
  if (file.exists(x$hdf5.fn)) {
    stop("HDF5 file already exists")
  }

  ## drop and initialize hdf5_xref table
  .init_hdf_xref_table(x)

  ## update gene_info table
  .update_gene_info_table(x)

  ## dump counts into hdf5 file and table
  .create_hdf5_expression_matrix(x)

  ## drop expression tables
  if (drop.expr.table) {
    dbRemoveTable(x$con, 'expression')
    dbGetQuery(x$con, 'DROP INDEX IF EXISTS expression_gene')
    dbGetQuery(x$con, 'DROP INDEX IF EXISTS expressoin_sample_gene')
    dbGetQuery(x$con, 'VACUUM')
  }
}

.create_hdf5_expression_matrix <- function(x) {
  stopifnot(is.FacileDataSet(x))
  datasets <- sample_stats_tbl(x) %>% collect(n=Inf)
  d.stats <- datasets %>%
    group_by(dataset) %>%
    summarize(n=n()) %>%
    ungroup

  cache <- file.path(x$parent.dir, 'count-cache')
  if (!dir.exists(cache)) {
    dir.create(cache)
    for (ds in d.stats$dataset) {
      message("===== caching ", ds, " =====")
      samples <- filter(datasets, dataset == ds)
      y <- fetch_expression.db(x, samples) %>% as.DGEList
      stopifnot(all.equal(y$genes$hdf5_index, 1:nrow(y)))
      saveRDS(y, file.path(cache, paste0(ds, '-DGEList.rds')))
    }
  }

  library(rhdf5)
  h5createFile(x$hdf5.fn)
  h5createGroup(x$hdf5.fn, 'expression')
  h5createGroup(x$hdf5.fn, 'expression/rnaseq')

  stuff <- lapply(d.stats$dataset, function(ds) {
    y <- readRDS(file.path(cache, paste0(ds, '-DGEList.rds')))
    cnts <- y$counts
    colnames(cnts) <- sub('^.*?_', '', colnames(cnts))
    name <- paste0('expression/rnaseq/', ds)
    message("==== creating ", name, " ====")
    h5createDataset(x$hdf5.fn, name, dim(cnts), storage.mode = "integer",
                    chunk=c(min(5000, nrow(cnts)), ncol(cnts)), level=4)
    h5write(cnts, file=x$hdf5.fn, name=name)
    tibble(dataset=ds, sample_id=colnames(cnts), hd5_index=1:ncol(cnts))
  })

  xref <- bind_rows(stuff)
  dbWriteTable(x$con, 'hdf5_sample_xref', xref, append=TRUE)
  xref
}

.init_hdf_xref_table <- function(x) {
  stopifnot(is.FacileDataSet(x))

  tbl.name <- 'hdf5_sample_xref'
  ## check if hdf5 xref tables exist
  qry <- paste0("SELECT name FROM sqlite_master WHERE type='table' AND ",
                "name='", tbl.name, "';")
  res <- dbGetQuery(x$con, qry)
  if (nrow(res) == 0) {
    ## table doesn't exist
  } else {
    dbRemoveTable(x$con, tbl.name)
  }

  qry <- paste(
    'CREATE TABLE', tbl.name,
    '(dataset TEXT, sample_id TEXT, hdf5_index INTEGER,',
    'PRIMARY KEY (dataset, sample_id));')

  dbGetQuery(x$con, qry)
}

.update_gene_info_table <- function(x) {
  ## add an hdf5_index column to gene_info table if it's not already there
  gi <- gene_info_tbl(x) %>%
    collect(n=Inf) %>%
    mutate(., hdf5_index=1:nrow(.))
  if (FALSE) {
    saveRDS(gi, file.path(x$parent.dir, 'gene-info.rds'))
  }

  ## Drop the table
  dbRemoveTable(x$con, 'gene_info')

  ## Re enter table with column
  table.sql <- paste0(
    "CREATE TABLE gene_info (",
    "feature_id TEXT, feature_type TEXT, symbol TEXT, n_exons INTEGER, ",
    "length INTEGER, source TEXT, hdf5_index INTEGER, ",
    "PRIMARY KEY (feature_id));")
  dbGetQuery(x$con, table.sql)
  dbGetQuery(x$con, 'DROP INDEX IF EXISTS gene_symbol')
  dbGetQuery(x$con, 'DROP INDEX IF EXISTS gene_feature_type')

  dbWriteTable(x$con, 'gene_info', gi, append=TRUE)
  dbSendQuery(x$con, 'CREATE INDEX gene_symbol ON gene_info (symbol);')
  dbSendQuery(x$con, 'CREATE INDEX gene_feature_type ON gene_info (feature_type);')
}

if (FALSE) {
  library(rhdf5)
  library(rprojroot)
  devtools::load_all(find_root(is_r_package))

  x <- FacileDataSet('~/workspace/data/facile/FacileDataSets/FacileAtezo-wip')
  createExpressionHDF5(x)
}

if (FALSE) {
  library(rhdf5)
  library(rprojroot)
  devtools::load_all(find_root(is_r_package))
  x <- FacileDataSet('~/workspace/data/facile/FacileDataSets/FacileTCGA-wip')
  createExpressionHDF5(x)
}


if (FALSE) {
  library(rhdf5)
  x <- FacileDataSet('/Users/lianogls/workspace/data/facile/FacileDataSets/FacileAtezo-wip')
  cache.dir <- file.path(x$parent.dir, 'count-cache')

  datasets <- sample_stats_tbl(x) %>% collect(n=Inf)
  d.stats <- datasets %>%
    group_by(dataset) %>%
    summarize(n=n()) %>%
    ungroup

  times <- lapply(d.stats$dataset, function(ds) {
    sy <- system.time(y <- readRDS(file.path(cache.dir, paste0(ds, '-DGEList.rds'))))
    s5 <- system.time(cnts <- h5read(x$hdf5.fn, paste0('expression/rnaseq/', ds)))
    yc <- y$counts
    dimnames(yc) <- NULL
    eq <- all.equal(cnts, yc)
    tibble(dataset=ds, rds.time=sy['elapsed'], hdf5.time=s5['elapsed'], all.equal=eq)
  })
  times <- bind_rows(times)
}

if (FALSE) {
  ## Test count equivalency -- seems to work
  fds <- FacileDataSet('~/workspace/data/facile/FacileDataSets/FacileAtezo-wip')
  genes <- c("800", "1009", "1289", "50509", "2191", "2335", "5159")
  samples <- sample_covariate_tbl(fds) %>%
    filter(variable == "indication" && value == 'bladder') %>%
    collect(n=Inf) %>%
    distinct(dataset, sample_id) %>%
    set_fds(fds)

  e.hdf5 <- fetch_expression(fds, samples, genes) %>%
    arrange(dataset, sample_id, feature_id) %>%
    collect
  e.sqlite <- fetch_expression.db(fds, samples, genes) %>%
    arrange(dataset, sample_id, feature_id) %>%
    collect
  all.equal(e.hdf5, e.sqlite) ## TRUE

  s2 <- samples %>% sample_n(20)
  e2.hdf5 <- fetch_expression(fds, s2, genes) %>%
    arrange(dataset, sample_id, feature_id) %>%
    collect
  e2.sqlite <- fetch_expression.db(fds, s2, genes) %>%
    arrange(dataset, sample_id, feature_id) %>%
    collect
  all.equal(e.hdf5, e.sqlite) ## TRUE

  system.time(all.exprs <- fetch_expression(fds, samples, as.matrix=TRUE))
  system.time(all.y <- as.DGEList(all.exprs))
  system.time(all.y2 <- as.DGEList(samples))
}
