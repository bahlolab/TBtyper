
#' @export
#' @importFrom SeqArray seqGetData seqSetFilter
#' @importFrom SeqVarTools variantInfo
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter arrange
#' @importFrom magrittr "%>%"
get_allele_counts_gds <- function(gds,
                                  var_info = NULL,
                                  ref = c('h37rv', 'mrca'),
                                  verbose = FALSE,
                                  na_to_zero = TRUE) {

  if (is.null(var_info)) {
    ref <- match.arg(ref)
    var_info <- get_var_info(ref = ref)
  }

  # check_args
  stopifnot(is_gds(gds),
            is.data.frame(var_info),
            setequal(c('variant_id', 'chr', 'pos', 'ref', 'alt'), colnames(var_info)))

  # find matching sites in gds
  var_id <- seqGetData(gds, 'variant.id')
  sam_id <- seqGetData(gds, 'sample.id')
  gr <- with(var_info, GenomicRanges::GRanges(chr, IRanges::IRanges(start = pos, width = 1L)))
  seqSetFilter(gds, gr, verbose = verbose)

  var_index <-
    variantInfo(gds, expand = TRUE) %>%
    as_tibble() %>%
    left_join(tibble(variant.id = seqGetData(gds, 'variant.id'),
                     num_allele = seqNumAllele(gds)),
              'variant.id') %>%
    inner_join(dplyr::rename(var_info, target = alt), by = c('chr', 'pos', 'ref')) %>%
    (function(x) {
      bind_rows(
        filter(x, allele.index == 1) %>% mutate(allele.index = 0) %>% dplyr::rename(allele = ref) %>% select(-alt),
        dplyr::rename(x, allele = alt) %>% select(-ref)
      )
    }) %>%
    arrange(variant.id, allele.index) %>%
    filter(allele.index == 0 | allele == target) %>%
    mutate(allele = if_else(allele.index == 0, 'ref', 'alt')) %>%
    select(-target) %>%
    spread(key = allele, value = allele.index) %>%
    arrange(variant.id) %>%
    mutate(ref.index = 1 + cumsum(num_allele) - num_allele,
           alt.index = ref.index + alt)

  # extract Allele Counts
  seqSetFilter(gds, variant.id = unique(var_index$variant.id), sample.id = sam_id, verbose = verbose)
  AD <- gds_get_AD_parallel(gds, verbose = verbose)$data

  ref_ac <- AD[, var_index$ref.index, drop = FALSE]
  alt_ac <- AD[, var_index$alt.index, drop = FALSE]

  if (na_to_zero) {
    alt_ac[is.na(alt_ac) & !is.na(ref_ac)] <- 0L
  }

  allele_counts <-
    array(c(ref_ac, alt_ac),
          dim = c(length(sam_id), nrow(var_index), 2),
          dimnames = list(sample = sam_id,
                          variant = var_index$variant_id,
                          allele = c('Ref', 'Alt')))

  return(allele_counts)
}


#' @importFrom magrittr "%>%"
#' @importFrom SeqArray seqGetData seqSetFilter
gds_get_AD_parallel <- function(gds, verbose = FALSE) {

  fn = gds$filename
  var.id = seqGetData(gds, 'variant.id')
  sam.id = seqGetData(gds, 'sample.id')
  workers = future::nbrOfWorkers()

  parallel::splitIndices(length(var.id), workers) %>%
    map( ~ var.id [.] ) %>%
    { .[lengths(.) > 0 ] } %>%
    furrr::future_map( ~{
      gds <- SeqArray::seqOpen(fn, allow.duplicate = T)
      seqSetFilter(gds, variant.id = ., sample.id = sam.id, verbose = verbose)
      seqGetData(gds, 'annotation/format/AD')
    }) %>%
    purrr::reduce(function(x, y) {
      list(length = c(x$length, y$length),
           data = cbind(x$data, y$data))
    })
}








