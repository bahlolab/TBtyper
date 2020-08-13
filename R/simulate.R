
# simulate allele counts from a mixture
# depth given by poisson distribution
simulate_sample_counts <- function(mix_gts, mix_prop,
                                   mean_depth = 100L,
                                   err_rate = 1/4,
                                   err_p = 0.01,
                                   as_data_df = TRUE,
                                   pseudocounts = FALSE) {
  if (FALSE) {
    .geno <- t(get_geno())
    nodes <- sample(ncol(.geno), 3)
    mix_gts <- .geno[, nodes]
    mix_prop <- sample(1000, 3) %>% { . / sum(.) }
    rm(.geno)
  }


  stopifnot(is.matrix(mix_gts),
            all(mix_gts == 1 | mix_gts == 0),
            is_proportion(mix_prop),
            ncol(mix_gts) == length(mix_prop),
            is_scalar_double(err_rate),
            is_scalar_double(err_p))


  result <-
    tibble(prob = colSums(t(mix_gts) * mix_prop),
           dp = rpois(length(prob), mean_depth)) %>%
    mutate(prob = if_else(rbinom(n(), 1, err_rate) == 1, prob * (1 - err_p) + (1-prob) * err_p, prob),
           bac = rbinom(n(), dp, prob),
           aac = dp - bac) %>%
    select(aac, bac, dp)

  if (as_data_df) {
    if (pseudocounts) {
      return(select(result, bac, dp) %>% mutate(bac = bac + 1, dp = dp + 2, err = 1 / dp))
    } else {
      return(select(result, bac, dp) %>% mutate(err = 0 ))
    }
  } else {
    return(cbind(result$aac, result$bac))
  }
}

#' @importFrom tidyr unnest
test_optim_phy_mix <- function(mix_n = 3L,
                               n_test = 100L,
                               mean_depth = 100L,
                               err_rate = 1/4,
                               err_p = 0.01,
                               pseudo = c('true', 'false', 'compare')) {

  pseudo <- match.arg(pseudo)
  geno <- t(get_geno())

  test_df <-
    tibble(test = seq_len(n_test)) %>%
    mutate(
      mix_def = map(test, ~ {
        tibble(mix_nodes = sample(ncol(geno), mix_n),
               mix_index = seq_along(mix_nodes),
               mix_prop = sample(1000, mix_n) %>% { . / sum(.) } )
      }),
      mix_data = map(mix_def, function(md) {
        simulate_sample_counts(mix_gts = geno[, md$mix_nodes],
                               mix_prop = md$mix_prop,
                               mean_depth = mean_depth,
                               err_rate = err_rate,
                               err_p = err_p,
                               pseudocounts = pseudo == 'true')
      }))

  result <-
    test_df %>%
    mutate(
      optim_mix = map2(mix_def, mix_data, function(mix_def, mix_data) {
        optim_phy_mix(data = mix_data, mix_gts = geno[, mix_def$mix_nodes]) %>%
          as_tibble() %>%
          select(est_prop = prop, n_iter)
      })
    ) %>%
    select(-mix_data) %>%
    unnest(c(mix_def, optim_mix))

  if (pseudo == 'compare') {
    result_2 <-
      test_df %>%
      mutate(
        optim_mix = map2(mix_def, mix_data, function(mix_def, mix_data) {
          optim_phy_mix(data = mutate(mix_data,
                                      bac = bac + 1,
                                      dp = dp + 2,
                                      err = 1 / dp),
                        mix_gts = geno[, mix_def$mix_nodes]) %>%
            as_tibble() %>%
            select(est_prop = prop, n_iter)
        })
      ) %>%
      select(-mix_data) %>%
      unnest(c(mix_def, optim_mix))

    result <- bind_rows(mutate(result, pseudocount = FALSE),
                        mutate(result_2, pseudocount = TRUE))
  }

  return(result)
}
