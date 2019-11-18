
is_proportion <- function(x, na.rm = T) {
  all((0 <= x) & (1 >= x), na.rm = na.rm)
}

#' @importFrom rlang is_scalar_double is_scalar_integer
is_scalar_proportion <- function(x, na.rm=T) {
  (is_scalar_double(x) || is_scalar_integer(x)) && is_proportion(x, na.rm = na.rm)
}


is_phylo <- function(x) {

  inherits(x, "phylo")
}

is_gds <- function(x) {

  inherits(x, "SeqVarGDSClass")
}
