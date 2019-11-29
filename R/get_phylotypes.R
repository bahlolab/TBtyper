
#' @export
#' @importFrom stringr str_c
get_geno <- function(ref = c('h37rv', 'mrca')) {

  ref <- match.arg(ref)
  readRDS(system.file(file.path(str_c(ref, '_phylotypes'), 'geno.rds'),
                      package = 'TBtyper'))
}

#' @export
#' @importFrom stringr str_c
get_phylo <- function(ref = c('h37rv', 'mrca')) {

  ref <- match.arg(ref)
  readRDS(system.file(file.path(str_c(ref, '_phylotypes'), 'phylo.rds'),
                      package = 'TBtyper'))
}

#' @export
#' @importFrom stringr str_c
get_var_info <- function(ref = c('h37rv', 'mrca')) {

  ref <- match.arg(ref)
  readRDS(system.file(file.path(str_c(ref, '_phylotypes'), 'var_info.rds'),
                      package = 'TBtyper'))
}
