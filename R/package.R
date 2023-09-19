#' @import checkmate
#' @import dplyr
#' @import methods
#' @importFrom utils read.csv
"_PACKAGE"

#' @importFrom broom tidy
#' @export
broom::tidy

# Export oft-used dplyr stuff --------------------------------------------------
# Should we just but dplyr in Depends?

#' @export
dplyr::arrange

#' @export
dplyr::collect

#' @export
dplyr::distinct

#' @export
dplyr::filter

#' @export
dplyr::group_by

#' @export
dplyr::mutate

#' @export
dplyr::rename

#' @export
dplyr::select

#' @export
dplyr::summarise

#' @export
dplyr::summarize

#' @export
dplyr::transmute

#' @export
dplyr::ungroup

#' @export
dplyr::left_join

#' @export
dplyr::inner_join

#' @export
dplyr::semi_join

#' @export
dplyr::anti_join

#' @export
dplyr::nest_join

#' @export
dplyr::bind_rows

#' @export
dplyr::bind_cols
