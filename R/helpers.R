
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


#' @importFrom rlang is_scalar_double is_scalar_integer
binom_likelihood <- function(x, size, p, err_01=0.005, err_10 = err_01, by_site = FALSE) {
  if (! is_proportion(p)) {
    warning('values outside [0,1] detected\n')
    p[p < 0] <- 0
    p[p > 1] <- 1
  }
  stopifnot(
    is_integerish(x),
    is_integerish(size),
    length(x) == length(size),
    length(x) == length(p) || length(p) == 1L,
    is_proportion(p),
    is_scalar_proportion(err_01),
    is_scalar_proportion(err_10)
  )
  p_ <- p * (1 - err_10) + (1 - p) * (err_01)
  res <- dbinom(x, size, p_, log = TRUE)
  if (by_site) {
    return(res)
  } else {
    return(sum(res, na.rm=TRUE))
  }
}

#' @importFrom dplyr bind_cols filter mutate select starts_with
#' @importFrom tidyr gather spread
#' @importFrom rlang is_double is_integer
fit_mix_prop <- function(data,
                         err_01 = 0.005,
                         err_10 = 0.005,
                         correct_baf = TRUE) {

  if (FALSE) {
    method = 'baf_est'
    data <- read_rds('test_data/phylomatch_data.rds')
    data$gts <- data$gts[,c(100, 200, 300)]
  }

  stopifnot(all(c('baf', 'gts') %in% colnames(data)),
            is_double(data$baf),
            is_integer(data$gts))

  nhap <- ncol(data$gts)

  if (correct_baf) {
    data$baf %<>% { (err_10 - .) / (err_01 + err_10 - 1) }
  }

  data %>%
    filter(!rowSums(gts) %in% c(0, nhap))  %>%
    { bind_cols(select(., -gts), as_tibble(.$gts %>% magrittr::set_colnames(str_c('H', seq_len(ncol(.)))))) } %>%
    mutate(site = 1:n()) %>%
    gather(starts_with('H'), key = 'ht', value = 'gt') %>%
    mutate(delta = gt - baf) %>%
    select(site, ht, delta) %>%
    spread(ht, delta) %>%
    select(-site) %>%
    as.matrix() %>%
    fit_proportions_QP()

}

fit_proportions_QP <- function(delta,
                               check_ident = F){

  if (check_ident) {
    d <- dist(t(delta)) %>% as.matrix()
    ident <- which(d == 0 & lower.tri(d), arr.ind = T)[,1]
    if (length(ident) > 0) {
      delta <- delta[, -ident, drop = F]
    }
  }
  n <- ncol(delta)

  chol_ <- tryCatch(
    chol(t(delta) %*% delta),
    error = function(e) {
      set.seed(0L)
      delta <- jitter(delta)
      chol(t(delta) %*% delta)
    })

  Rinv <- solve(chol_)
  C <- cbind(rep(1, n), diag(n))
  b <- c(1, rep(0, n))
  d <- rep(0, n)
  fit <- quadprog::solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
  # numerical rounding can result in values slightly outside contraints
  # correct values outside counstrains and ensure sum to 1
  sol <- fit$solution
  sol[sol < 0] <- 0
  sol %<>% { . / sum(.) }
  sol[sol > 1] <- 1

  if (check_ident) {
    if (length(ident) > 0) {
      sol_ <- rep(0, n + length(ident))
      sol_[-ident] <- sol
      sol <- sol_
    }
  }
  return(sol)
}
