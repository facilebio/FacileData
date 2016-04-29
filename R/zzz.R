.onLoad <- function(libname, pkgname) {
  opts <- options()

  ## This package serves as an "abstract implementation" to a FacileDb database.
  ## The packages that implement acces to a FacileWareshouse should define
  ## the following options:
  ##
  ##   - *.datapath
  ##   - *.dbpath
  ##   - *.cachedir

  ## The code below is the same in the FacileAtezo db package, because we will
  ## use that dataset for testing.
  if (Sys.getenv("HOME") == "/Users/lianogls") {
    dpath <- '/Users/lianogls/workspace/data/facile/atezo/v01'
  } else {
    dpath <- '/gne/devsci/data/obdroot/PDL1mab/pcd4989g/ngs/intermdata/lianoglou/facile-atezo/v01'
  }

  impl.prefix <- 'fatezo'
  dpath <- getOption(sprintf('%s.datapath', impl.prefix), dpath)

  pkg.opts <- list(
    datapath=dpath,
    dbpath=file.path(dpath, 'AtezoDb.sqlite'),
    cachedir=file.path(dpath, 'cache'),
    covdef=file.path(dpath, 'sample-meta-definitions.yaml'))
  names(pkg.opts) <- sprintf('%s.%s', impl.prefix, names(pkg.opts))
  toset <- !(names(pkg.opts) %in% names(opts))
  if (any(toset)) {
    options(pkg.opts[toset])
  }

  ## Check options
  db.path <- getOption(paste0(impl.prefix, '.dbpath'))
  if (!file.exists(db.path)) {
    msg <- paste0(
      "Default path to faciledb is not a valid file: ", db.path, "\n",
      "Set options('facile.datapath') before loading the facilewarehouse ",
      "package to a valid path to the SQLite database to skip this message.\n",
      "A good place to do this for your local work is in your ~/.Rprofile")
    warning(msg, immediate.=TRUE)
  }

  invisible()
}
