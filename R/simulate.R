
# simulate allele counts from a mixture
# depth given by poisson distribution
simulate_sample_counts <- function(mix_gts, mix_prop,
                                   mean_depth = 100L,
                                   err_rate = 1/4,
                                   err_p = 0.01,
                                   as_data_df = TRUE) {
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
    return(select(result, bac, dp) %>% mutate(bac = bac + 1, dp = dp + 2, err = 1 / dp))
  } else {
    return(cbind(result$aac, result$bac))
  }
}

test_optim_phy_mix <- function(mix_n = 3L,
                               n_test = 100L,
                               mean_depth = 100L,
                               err_rate = 1/4,
                               err_p = 0.01) {
  geno <- t(get_geno())

  test_df <-
    tibble(test = seq_len(n_test)) %>%
    mutate(
      mix_def = map(test, ~ {
        tibble(mix_nodes = sample(ncol(geno), mix_n),
               mix_index = seq_along(mix_nodes),
               mix_prop = sample(1000, mix_n) %>% { . / sum(.) } %>% sort())
      }),
      mix_data = furrr::future_map(mix_def, function(md) {
        simulate_sample_counts(mix_gts = geno[, md$mix_nodes],
                               mix_prop = md$mix_prop,
                               mean_depth = mean_depth,
                               err_rate = err_rate,
                               err_p = err_p)
      }),
      optim_mix = furrr::future_map2(mix_def, mix_data, function(mix_def, mix_data) {
        optim_phy_mix(data = mix_data, mix_gts = geno[, mix_def$mix_nodes]) %>%
          as_tibble() %>%
          select(est_prop = prop, n_iter)
      })
    ) %>%
    select(-mix_data) %>%
    unnest(c(mix_def, optim_mix))

  return(test_df)

  # test_df %>%
  #   select(-mix_data) %>%
  #   unnest(c(mix_def, optim_mix)) %>%
  #   ggplot(aes(mix_prop, est_prop)) +
  #   geom_point(alpha = 0.25) +
  #   facet_wrap(~mix_index)
}
