
#' @export
get_geno <- function() {

  readRDS(system.file(file.path('h37rv_phylotypes', 'geno.rds'),
                      package = 'TBtyper'))
}

#' @export
get_phylo <- function() {

  readRDS(system.file(file.path('h37rv_phylotypes', 'phylo.rds'),
                      package = 'TBtyper'))
}

#' @export
get_var_info <- function() {

  readRDS(system.file(file.path('h37rv_phylotypes', 'var_info.rds'),
                      package = 'TBtyper'))
}

