
#' @export
#' @importFrom SeqArray seqGetData seqSetFilter
#' @importFrom SeqVarTools variantInfo
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter
#' @importFrom magrittr "%>%"
get_allele_counts_gds <- function(gds,
                                  var_info = NULL,
                                  ref = c('h37rv', 'mrca'),
                                  verbose = FALSE,
                                  alt_na_to_zero = TRUE) {

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

  var_match <-
    variantInfo(gds, expanded = TRUE) %>%
    as_tibble() %>%
    group_by(variant.id) %>%
    mutate(num_alt = max(allele.index)) %>%
    ungroup() %>%
    inner_join(var_info, c('chr', 'pos', 'ref', 'alt')) %>%
    filter(variant.id %in% var_id) %>%
    distinct()

  # extract Allele Counts
  seqSetFilter(gds, variant.id = var_match$variant.id, sample.id = sam_id, verbose = verbose)
  AD <- gds_get_AD_parallel(gds, verbose = verbose)$data
  ri <- with(var_match, cumsum(num_alt + 1) - (num_alt))
  ai <- (ri + var_match$allele.index) %>% na.omit() %>% c()

  RAC <- AD[, ri, drop = FALSE]
  AAC <- AD[, ai, drop = FALSE]

  if (alt_na_to_zero) {
    AAC[is.na(AAC) & !is.na(RAC)] <- 0L
  }


  allele_counts <-
    array(c(RAC, AAC),
          dim = c(length(sam_id), length(ri), 2),
          dimnames = list(sample = sam_id,
                          variant = var_match$variant_id,
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








