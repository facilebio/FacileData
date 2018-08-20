# Cleanup Branch

This includes changes that:

1. Lays some groundwork down for factoring out the "FcileData API" from the
   FacileDataSet.
2. Largely fills out documentation needed to placate `R CMD check` [not done]
3. Fills out some vignettes [not done]
4. pkgdown

## Introduction of FacileDataStore

Introduced a `FacileDataStore` "abstract class" in anticipation of refactoring
out the "FacileData API" into a top-level `FacileData` package.

The idea is that:
  
  1. Any object that implements the "FacileData API" must include
     `"FacileDataStore"` in its class hierarchy (at the root(?)). For instance,
     `class(exampleFacileDataSet())` returns
     `"ExampleFacileTCGADataSet" "FacileDataSet" "FacileDataStore"`
  
  2. The S3 methods that will be factored out of this package to define the
     "FacileData API" will effectively use `*.FacileDataStore` as their
     "default" methods. `*.default` FacileData API S3 methods should either
     (i) not be defined; or (ii) throw an error.
  
Note that the names of the to-be-factored-out base package name ("FacileData",
here) and "FacileData API" are up for discussion. I'm just using them here
as placeholders to reference the concept we are all working towards.

Random notes in this orbit:

* I added `@family API` roxygen tags to methods that I (loosely) think should
  make up the "FacileData API". I'm pretty sure

* `@family API` methods should probably `assert_facile_data_store()` instead 
  of `assert_facile_data_set()`

## FacileDataSet validity checking

This is currently split across many functions:

  * `(assert|check|test)_facile_data_set()` are [checkmate][checkmate]-esque
  * The `FacileDataSet()` constructor calls `validate.facile.dirs()` which
    is a workhorse-of a function.
  * `is.FacileDataSet()` now delegates to `assert_facile_data_set()`

ASK: Should `check_facile_data_set()` do all of the checking that
`validate.facile.dirs()` does? The downside is that `assert_facile_data_set()`
is likely called a lot, and `validate.facile.dirs()` may be doing too much
all of the time?

[checkmate]: https://CRAN.R-project.org/package=checkmate

## Minor Changes

* Enables roxygen markdown parsing as default in DESCRIPTION.
* Most `##'` documentation blocks are changed to `#'`

