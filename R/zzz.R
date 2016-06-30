.onLoad <- function(libname, pkgname) {
  impl.prefix <- 'ftest'
  ## This package serves as an "abstract implementation" to a FacileDb database.
  ## The packages that implement acces to a FacileWareshouse should define
  ## the following options:
  ##
  ##   - *.datapath
  ##   - *.dbpath
  ##   - *.covdef
  ##   - *.cachedir
  ##
  ## Since this FacileRepo package should not be tied to a specific FacileRepo
  ## "implementation", we set the impl.prefex to be "ftest". This means that
  ## this setup will create (or reuse) the following global options
  ##
  ##   - ftest.datapath
  ##   - ftest.dbpath
  ##   - ftest.covdef
  ##   - ftest.cachedir

  ## Default database this will point to is configured so everything works
  ## when this package is loaded (deployed) on rescomp
  ## TODO: Update this database path to a test db instead of Atezo

  ## Until we create a test database and distribute within the package, I'm
  ## affraid we can't avoid explicity defining the *.datapath. This is because
  ## the unit tests are run in a "clean" (R --vanilla) environment which doesn't
  ## load the stuff in your .Rprofile
  dpath <- system.file('extdata', 'test', package='FacileRepo')
  db.name <- 'TcgaDb-test.sqlite'
  dpath <- getOption(sprintf('%s.datapath', impl.prefix), dpath)

  pkg.opts <- list(
    datapath=dpath,
    dbpath=file.path(dpath, db.name),
    cachedir=file.path(dpath, 'cache'),
    covdef=file.path(dpath, 'sample-meta-definitions.yaml'))
  names(pkg.opts) <- sprintf('%s.%s', impl.prefix, names(pkg.opts))

  ## We only set these options if they aren't already set in the global options
  ## The developers should set the appropriate options in the ~/.Rprofile
  opts <- options()
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
    ## warning(msg, immediate.=TRUE)
  }

  invisible()
}
