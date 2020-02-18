.biocboxes <- c("DGEList", "EList", "ExpressionSet")

#' Assembles a Bioconductor data container for a given assay.
#'
#' There is a default bioc class provided for different assay types, however
#' the class type can be overrided by the `class` parameter.
#'
#' TODO: Migrate `as.DGEList` functionality to `biocbox()`
#' @noRd
#' @export
biocbox.facile_frame <- function(x, assay_name = default_assay(x),
                                 features = NULL, class = NULL,
                                 custom_key = Sys.getenv("USER"), ...) {
  class <- match.arg(class)
  biocbox(fds(x), assay_name = assay_name, )
}

biocbox.FacileDataStore <- function(x, assay_name = default_assay(x),
                                    features = NULL,
                                    samples = NULL, class = NULL,
                                    custom_key = Sys.getenv("USER"), ...) {
  class <- match.arg(class)
}
