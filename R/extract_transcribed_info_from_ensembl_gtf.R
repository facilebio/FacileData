# This file has functions to call over a genocde/ensembl GTF to get
# transcript- and gene-level meta information for the stuff
# It doesn't really belong in FacileData proper, but is useful when generating
# new datasets "at large"

#' Extract gene- and transcript-level information from an ENSEMBL gtf.
#'
#' This was written for release_28 annotations. This is noted because some
#' column names seemsed to have changed, ie. "gene_type" instead of
#' "gene_biotype", etc. Let's see how consistent this is!
#'
#' @export
#'
#' @param fn the path to the ENSEMBL (or GENCODE) GTF
#' @return a list of tibbles with `$transcript_info` and `$gene_info` elements
extract_transcribed_info_from_ensembl_gtf <- function(
    fn, gene_type = "gene_type", transcript_type = "transcript_type") {
  assert_string(gene_type)
  assert_string(transcript_type)

  if (!requireNamespace("rtracklayer")) {
    stop("rtacklayer required to play with gtfs")
  }
  if (!requireNamespace("S4Vectors")) {
    stop("S4Vectors required to play with gtfs")
  }

  gtf <- rtracklayer::import(fn)

  # Working with at least release_28?
  assert_choice(gene_type, colnames(S4Vectors::mcols(gtf)))
  S4Vectors::mcols(gtf)[[gene_type]] <-
    .level_biotypes(S4Vectors::mcols(gtf)[[gene_type]])

  if (transcript_type != gene_type) {
    assert_choice(transcript_type, colnames(S4Vectors::mcols(gtf)))
    S4Vectors::mcols(gtf)[[transcript_type]] <-
      .level_biotypes(S4Vectors::mcols(gtf)[[transcript_type]])
  }

  tx.info <- .transcript_info(gtf, transcript_type)
  gn.info <- .gene_info(gtf, tx.info, gene_type)

  # Note that stripping the version information from gene_id (and transcript_id)
  # will result in duplicate identifiers. These identifiers come from the
  # pseudoautosomal regions of the X and Y chromosomes.
  list(transcript_info = tx.info, gene_info = gn.info)
}

#' Utility function to "factorize" biotypes into an order we care about.
#'
#' ENSEMBL GTFs provide biotype information for genes/transcripts. These are
#' things like "3prime_overlapping_ncRNA", "antisense", ..., "protein_coding",
#' etc. This function turns the "biotype"-vector `x` into a factor with levels
#' in (roughly) the order we care to "unique"-ify these levels. Ie. if a gene
#' has a "protein_coding" annotation, we will care to keep that one over one
#' of its annotations which categorize it as a "processed_transcript"
#'
#' @param x a `character` vector of biotypes
#' @return a factor version of `x`, with `levels(x)` in approximately the order
#'   we care about.
.level_biotypes <- function(x) {
  assert_character(x)
  # Set the order of gene/transcript "type" (formerly known as biotype) to be
  # approximately what we find most informative when summarizing multiple
  # transcripts to the gene level.
  # types <- unique(c(unique(gtf$gene_type), unique(gtf$transcript_type)))
  types <- unique(x)
  bt.order <- c('protein_coding', 'lincRNA',
                types[grepl("lncrna", types, ignore.case = TRUE)],
                'miRNA', 'rRNA', 'snoRNA',
                'scRNA', 'scaRNA', 'sRNA' ,"antisense", "ERCC")
  bt.order <- c(bt.order, setdiff(types, bt.order))
  bt.order <- intersect(bt.order, types)
  out <- factor(x, bt.order)
  if (any(!is.na(x) & is.na(out))) {
    stop(".level_biotypes messed up somehow")
  }
  out
}

#' Helper function that creats augmented transcript information file.
#'
#' @noRd
#'
#' @param gr a GRanges object initialized from an ENSEMBL/GENOCDE GTF
#' @return a table of transcript information
.transcript_info <- function(gr, transcript_type) {
  requireNamespace("GenomicRanges", quietly = TRUE)
  dfa <- GenomicRanges::as.data.frame(gr)

  tx.length <- dfa %>%
    filter(type == "exon") %>%
    group_by(transcript_id) %>%
    summarize(length = sum(width))

  tx.meta <- dfa %>%
    filter(type == "transcript") %>%
    select(transcript_id, transcript_name, {{transcript_type}},
           gene_name, gene_id, source, seqnames, start, end, strand)
  if (any(duplicated(tx.meta$transcript_id))) {
    stop("Duplicated transcript ids found in gtf file")
  }
  if (!setequal(tx.length$transcript_id, tx.meta$transcript_id)) {
    stop("tx.meta and tx.length are not setequal")
  }

  inner_join(tx.meta, tx.length, by = "transcript_id")
}

#' Helper function to create augmented gene-level metadata.
#'
#' @noRd
#'
#' @param gr a GRanges object initialized from an ENSEMBL/GENOCDE GTF
#' @param tx_info the tibble that came from .transcript_info()
#' @return a tibble of gene-level information
.gene_info <- function(gr, tx_info, gene_type) {
  requireNamespace("GenomicRanges", quietly = TRUE)
  gr.exons <- gr[gr$type == "exon"]

  grl <- GenomicRanges::split(gr.exons, gr.exons$gene_id)
  rgrl <- GenomicRanges::reduce(grl)
  gwidths <- sum(GenomicRanges::width(rgrl))
  gwidths <- tibble(gene_id = names(gwidths), length = unname(gwidths))

  dfa <- as.data.frame(gr)

  # Arrange genes by gene_id,gene_type,seqnames to take the "most useful guess"
  # for its metadata
  g.info <- dfa %>%
    filter(type == "gene") %>%
    arrange(gene_id, source, gene_type, seqnames) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    select(gene_id, symbol = gene_name, {{gene_type}}, seqnames, start, end,
           strand, source)

  ntx <- tx_info %>%
    group_by(gene_id) %>%
    summarize(ntxs = n())

  out <- g.info %>%
    left_join(gwidths, by = "gene_id") %>%
    left_join(ntx, by = "gene_id") %>%
    select(gene_id:gene_type, ntxs, length, everything())

  out
}
