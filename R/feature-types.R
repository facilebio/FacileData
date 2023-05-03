#' Guesses the type of feature identifiers from a character vector.
#'
#' We rely on meta-information about our data types than "usual", and its useful
#' to know what types of identifiers we are using for different assay. This
#' function tries to guess whether an identifier is an ensembl gene identifier,
#' entrez id, etc.
#'
#' A two-column data.frame is returned for id_type and organism. Organism
#' is "unknown" for identifiers where there this can't be inferred (like Refseq).
#'
#' If an identifier matches more than one id_type, the id_type is set to
#' `"ambiguous"`. If the identifier doesn't match any guesses, then `"unknown"`.
#'
#' @export
#' @param x a character vector of ids
#' @return data.frame with `id` (`x`) and `id_type`. If `with_organism = TRUE`,
#'   a third `organism` column is added with a guess for the organism.
#' @examples
#' fids <- c("NC_000023", "ENSG00000101811", "ENSMUSG00000030088.2", "85007")
#' infer_feature_type(fids)
infer_feature_type <- function(x, with_organism = FALSE, ...) {
  regex <- list(
    # ens_gene = "^ENS[A-Z]*G\\d+(\\.\\d+)?$",
    # ens_tx   = "^ENS[A-Z]*?T\\d+(\\.\\d+)?$",
    ensgid   = "^ENS[A-Z]*G\\d+(\\.\\d+)?$",
    enstid   = "^ENS[A-Z]*?T\\d+(\\.\\d+)?$",
    refseq   = "^[NXW][CGMRP]_\\d+(\\.\\d+)?$",
    entrez   = "^\\d+$")

  bool <- sapply(regex, grepl, x)
  nmatch <- rowSums(bool)
  type <- names(regex)[apply(bool, 1, function(vals) which(vals)[1])]
  type <- ifelse(nmatch == 1L, type, "ambiguous")
  type <- ifelse(nmatch == 0L, "unknown", type)

  is.bad <- type %in% c("ambiguous", "unknown")
  if (any(is.bad)) {
    warning(sum(is.bad), " identifiers were either ambiguous or unknown",
            immediate. = TRUE)
  }

  out <- tibble(
    id = x,
    id_type = type)

  if (with_organism) {
    is.ens <- grepl("^ens_", out[["id_type"]])
    ens <- sub("^ENS", "", out[["id"]])
    is.human <- is.ens & grepl("^[TG]\\d+", ens)
    is.mouse <- is.ens & grepl("MUS[TG]\\d+", ens)
    out[["source_organism"]] <- ifelse(is.human, "Homo sapiens", "unknown")
    out[["source_organism"]] <- ifelse(is.mouse, "Mus musculus", out[["source_organism"]])
    unk <- out[["source_organism"]] == "unknown"
    if (any(unk)) {
      warning(sum(unk), " identifiers could not be matched to an organism",
              immediate. = TRUE)
    }
  }

  out
}
