

# allele_counts should be a three dimensional array
# rows == sample, cols == variant, dim 3-1 == ref_allele count, dim 3-2 == alt_allele count
# e.g allele_counts <- array (c(ref_count, alt_count), dim c (dim(ref_count), 2))
#' @export
#' @importFrom rlang is_integer is_scalar_integerish is_bool
#' @importFrom treeio Nnode2
#' @importFrom magrittr set_colnames set_rownames
fit_phylotypes <- function(allele_counts,
                           phylo = NULL,
                           geno = NULL,
                           ref = c('h37rv', 'mrca'),
                           max_phylotypes = 5L,
                           min_mix_prop = 0.005,
                           min_depth = 50L,
                           max_depth = 250L,
                           max_p_val_perm = 0.001,
                           max_p_val_lrt = 0.001,
                           spike_in_p = 0.01,
                           reoptimise = TRUE,
                           reoptimise_max_iter = 5L,
                           min_sites = 1000L,
                           n_perm = 1000L,
                           exclude_parent = TRUE,
                           exclude_child = TRUE,
                           exclude_ancestor = TRUE,
                           exclude_descendant = FALSE,
                           exclude_distance = 50L,
                           exclude_inner = FALSE
                           ) {
  ref <- match.arg(ref)

  if (is.null(phylo)) {
    phylo <- get_phylo(ref = ref)
  }

  if (is.null(geno)) {
    geno <- get_geno(ref = ref)
  }

  # check args
  stopifnot(
    is_integer(allele_counts) && is.array(allele_counts) && min(allele_counts, na.rm = T) >= 0L,
    length(dim(allele_counts)) == 3L && dim(allele_counts)[3] == 2L,
    is_integer(geno) && is.matrix(geno) && max(geno, na.rm = T) == 1L && min(geno, na.rm = T) == 0L,
    is_phylo(phylo),
    setequal(rownames(geno), c(phylo$tip.label, phylo$node.label)),
    all(colnames(allele_counts) %in% colnames(geno)),
    is_scalar_integerish(max_phylotypes) & max_phylotypes > 0,
    is_scalar_integerish(min_depth) & min_depth > 0,
    is_scalar_integerish(max_depth) & max_depth > 0,
    is_scalar_proportion(max_p_val_perm),
    is_scalar_proportion(max_p_val_lrt),
    is_scalar_integerish(n_perm))

  # subset genotypes to those in input allele counts
  geno_sub <- geno[, colnames(allele_counts)]
  message('Note: using ', ncol(geno_sub), ' of ', ncol(geno), ' genotypes')

  # subset phylotypes to those that can be distinguished over input genotypes
  phylo_sub <- collapse_phylotypes(phylo, geno_sub)
  message('Note: using ', Nnode2(phylo_sub), ' of ', Nnode2(phylo), ' phylotypes')

  # transpose genotypes
  t_geno <- t(geno_sub)[, node_to_label(phylo_sub, seq_len(Nnode2(phylo_sub)))]
  rate <- rowSums(t_geno) / ncol(t_geno)

  perm_geno <-
    vapply(seq_len(nrow(t_geno)),
           function (i) { rbinom(n_perm, 1, rate[i]) },
           integer(n_perm)) %>%
    t()

  # perm_geno <- matrix(0L, nrow = nrow(t_geno), ncol = ncol(t_geno))
  # for (i in seq_len(nrow(t_geno))) {
  #   perm_geno[i, sample(ncol(t_geno), sum(t_geno[i, ]))] <- 1L
  # }

    vapply(rowSums(t_geno),
           function(i) sample(ncol(t_geno), i),
           integer(ncol(t_geno)))

  results <-
    seq_len(nrow(allele_counts)) %>%
    furrr::future_map_dfr(function(i) {
      fit_sample(
        phylo = phylo_sub,
        gts = t_geno,
        pgts = perm_geno,
        sm_allele_counts = allele_counts[i, ,],
        max_phylotypes = max_phylotypes,
        min_mix_prop = min_mix_prop,
        min_depth = min_depth,
        max_depth = max_depth,
        max_p_val_perm = max_p_val_perm,
        max_p_val_lrt = max_p_val_lrt,
        spike_in_p = spike_in_p,
        reoptimise = reoptimise,
        reoptimise_max_iter = reoptimise_max_iter,
        min_sites = min_sites,
        exclude_parent = exclude_parent,
        exclude_child = exclude_child,
        exclude_ancestor = exclude_ancestor,
        exclude_descendant = exclude_descendant,
        exclude_distance = exclude_distance,
        exclude_inner = exclude_inner) %>%
        mutate(sample_id = rownames(allele_counts)[i])
    }) %>%
    select(sample_id, tidyr::everything())

  return(results)
}

#' @importFrom tidyr replace_na chop
#' @importFrom magrittr "%<>%"
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter if_else bind_rows tibble first last
fit_sample <- function(phylo,
                       gts,
                       pgts,
                       sm_allele_counts,
                       max_phylotypes,
                       min_mix_prop,
                       min_depth,
                       max_depth,
                       max_p_val_perm,
                       max_p_val_lrt,
                       spike_in_p,
                       reoptimise,
                       reoptimise_max_iter,
                       min_sites,
                       exclude_parent,
                       exclude_child,
                       exclude_ancestor,
                       exclude_descendant,
                       exclude_distance,
                       exclude_inner,
                       fuzzy_phylotypes = TRUE,
                       fuzzy_max_dist = 100L) {

  # TODO - bring back optim_phylotype as slide_phylotype (could be anywhere between current node, parents and children)
  #   - give optimised a name like 1.1.2/2_0.34_1.1.2/2/1

  # bac = b-allele count, dp = depth, gts = phylotype genotypes, pgts = permutation genotypes
  data <-
    tibble(variant = rownames(sm_allele_counts),
           bac = sm_allele_counts[,'Alt'],
           dp = rowSums(sm_allele_counts),
           gts = gts,
           pgts = pgts) %>%
    filter(dp >= min_depth) %>%
    # round down samples with depth greater than max_depth
    mutate(dp_gt_mx = dp > max_depth,
           bac = replace(bac, dp_gt_mx, round(max_depth * (bac / dp)[dp_gt_mx]) %>% as.integer()),
           dp = replace(dp, dp_gt_mx, max_depth)) %>%
    select(-dp_gt_mx) %>%
    # add pseudocounts and calculate 'error' rate for binomial model
    mutate(bac = bac + 1L,
           dp = dp + 2L,
           err = 1 / dp) %>%
    na.omit()

  phy_gts <- data$gts %>% set_colnames(node_to_label(phylo, seq_len(Nnode2(phylo)))) %>%
    set_rownames(data$variant)
  perm_gts <- data$pgts
  data %<>% select(-gts, -pgts)

  # arbitrary threshold on min sites
  if (nrow(data) < min_sites) {  return(tibble(note = "insufficient data"))  }

  node_dist <- phylo_geno_dist(phylo, t(phy_gts))

  # find the first phylotype
  res_1 <-
    score_phy_mix(data, phy_gts, perm_gts) %>%
    rename(node = index) %>%
    arrange(desc(likelihood), p_val_perm) %>%
    filter(! node %in% exclusions(integer(), phylo, node_dist,
                                  exclude_inner = exclude_inner)) %>%
    mutate(p_val_lrt = pchisq(-2 * (likelihood - likelihood[1]),
                              node_dist[node[1], node],
                              lower.tail = F) %>% replace(1, NA))

  res_1_top <-
    res_1 %>%
    filter(p_val_perm < max_p_val_perm) %>% # rho > min_rho) %>%
    slice(1)

  # store the best match in sample_fit
  if (nrow(res_1_top) > 0) {

    sample_fit <-
      res_1_top %>%
      mutate(mix_n = 1, mix_index = 1, mix_prop = 1, search = list(res_1),
             p_val_lrt = 0, fuzzy_sites = list(integer())) %>%
      select(mix_n, mix_index, mix_prop, node, phylotype, mix_prop, likelihood,
             p_val_perm, p_val_lrt, search, fuzzy_sites)

    if (fuzzy_phylotypes) {
      fuzz <- optim_fuzzy(data, phy_gts = phy_gts, phylo = phylo, node_dist = node_dist,
                          max_dist = fuzzy_max_dist, mix_nodes = sample_fit$node[1], mix_fit = 1)
      sample_fit[1,] %<>% mutate(likelihood = likelihood + fuzz$lh_delta,
                                 fuzzy_sites = list(fuzz$sites))
    }

  } else {
    return(tibble(note = "no matches passing filters"))
  }

  for (i in seq_len(max_phylotypes - 1)) {
    lh_last <- filter(sample_fit, mix_n == i) %>% pull(likelihood) %>% first()
    nodes_last <- filter(sample_fit, mix_n == i) %>% pull(node)
    mix_prop_last <- filter(sample_fit, mix_n == i) %>% pull(mix_prop)
    mix_gts_last <-
      filter(sample_fit, mix_n == i) %>%
      pmap(function(node, fuzzy_sites, ...) flip_at(phy_gts[, node, drop = FALSE], fuzzy_sites)) %>%
      do.call('cbind', .)

    # find next mixture component
    res_next <-
      score_phy_mix(data = data,
                    phy_gts_search = phy_gts[, -nodes_last, drop = FALSE],
                    perm_gts = perm_gts,
                    phy_gts_fixed = mix_gts_last,
                    fit = c_spike_in(mix_prop_last, spike_in_p)) %>%
      rename(node = index) %>%
      arrange(desc(likelihood), p_val_perm) %>%
      filter(!node %in% exclusions(nodes_last, phylo, node_dist,
                                   exclude_parent = exclude_parent,
                                   exclude_child = exclude_child,
                                   exclude_ancestor = exclude_ancestor,
                                   exclude_descendant = exclude_descendant,
                                   exclude_distance = exclude_distance,
                                   exclude_inner = exclude_inner))

    node_next <- res_next %>% filter(p_val_perm < max_p_val_perm) %>% pull(node) %>% first()
    # stop if we havent found a new candidate node
    if (length(node_next) != 1) { break }

    n_diff_next <-
      map_int(seq_len(ncol(mix_gts_last)), ~ sum(mix_gts_last[, .] != phy_gts[, node_next])) %>%
      sum(na.rm = TRUE)

    optim_prop <- optim_phy_mix(data, mix_gts = cbind(mix_gts_last, phy_gts[, node_next]))
    lrt_p_val <- pchisq(-2 * (lh_last - optim_prop_0$likelihood), n_diff_next, lower.tail = F)

    mix_fit <- with(optim_prop,
                    tibble(mix_n = i + 1,
                           mix_index = seq_len(mix_n),
                           mix_prop = prop,
                           node = c(nodes_last, node_next),
                           search = map(node, ~ tibble(phylotype = character(),
                                                       p_val_perm = numeric())),
                           likelihood = likelihood,
                           p_val_lrt = lrt_p_val))

    mix_fit$search[[length(mix_nodes)]] <- res_next

    # reoptimise existing mixture components in light of new component
    if (reoptimise && lrt_p_val < max_p_val_lrt) {
      n_iter <- 0L
      n_alt <- 0L
      n_alt_last <- rep(-1L, i+1)
      n_alt_last[i+1] <- 0L
      while(n_iter < reoptimise_max_iter) {
        n_iter <- n_iter + 1L
        change_made <- FALSE
        for (j in seq_along(mix_nodes)) {
          if (n_alt > n_alt_last[j]) {
            n_alt_last[j] <- n_alt
            # find highest scoring alternate components
            mix_fit$search[[j]] <-
              score_phy_mix(phylo = phylo,
                            data = data,
                            nodes_fixed = mix_fit$node[-j],
                            fit = c(mix_fit$mix_prop[-j], mix_fit$mix_prop[j])) %>%
              arrange(desc(likelihood), p_val_perm) %>%
              filter(!node %in% exclusions(mix_fit$node[-j], phylo, node_dist,
                                           exclude_parent = exclude_parent,
                                           exclude_child = exclude_child,
                                           exclude_ancestor = exclude_ancestor,
                                           exclude_descendant = exclude_descendant,
                                           exclude_distance = exclude_distance,
                                           exclude_inner = exclude_inner)) %>%
              arrange(desc(likelihood), p_val_perm) %>%
              mutate(lrt_p_val = pchisq(-2 * (likelihood - likelihood[1]),
                                        node_dist[node[1], node],
                                        lower.tail = F) %>% replace(1, NA))
              # mutate(lrt_p_val = pchisq(-2 * (likelihood - likelihood[1]), nrow(data), lower.tail = F))

            top <-
              mix_fit$search[[j]] %>%
              filter(p_val_perm < max_p_val_perm,
                     node != mix_fit$node[j]) %>%
              slice(1)

            if (nrow(top) > 0) {
              # find maximum likelihood proportion
              optim_prop <- optim_phy_mix(data,
                                          nodes = replace(mix_fit$node, j, top$node),
                                          fit = mix_fit$mix_prop)
              # accept is proportion is highed than curren
              if (optim_prop$likelihood > mix_fit$likelihood[1]) {
                mix_fit %<>%
                  mutate(likelihood = optim_prop$likelihood,
                         node = optim_prop$nodes,
                         mix_prop = optim_prop$prop)
                change_made <- TRUE
                n_alt <- n_alt + 1L
              }
            }
          }
        }
        # exit if we couldn't find a better solution
        if (!change_made) { break }
      }
    }

    mix_fit %<>%
      mutate(phylotype = node_to_label(phylo, node),
             p_val_perm = map_dbl(search, ~ .$p_val_perm[1]))

    sample_fit <- bind_rows(sample_fit, mix_fit)

    min_prop <- filter(sample_fit, mix_n == i+ 1) %>% pull(mix_prop) %>% min()
    if (min_prop < min_mix_prop || mix_fit$p_val_lrt[1] > max_p_val_lrt) { break }
  }

  return(chop(sample_fit, c(mix_prop, mix_index,node, phylotype, p_val_perm, search)))
}

#' @importFrom purrr map map_dbl pmap_dbl map2_dbl
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter summarise
score_phy_mix <- function(data, phy_gts_search, perm_gts,
                          phy_gts_fixed = matrix(0L, nrow(phy_gts_search), 0L),
                          fit = rep(1/(ncol(phy_gts_fixed) + 1), ncol(phy_gts_fixed) + 1)) {

  #check args
  stopifnot(is.data.frame(data),
            is.matrix(phy_gts_search),
            is.matrix(phy_gts_fixed),
            is.matrix(perm_gts),
            nrow(data) == nrow(phy_gts_search),
            nrow(data) == nrow(phy_gts_fixed),
            nrow(data) == nrow(perm_gts),
            is_double(fit) && all(fit >= 0 | fit <= 1) && (abs(1 - sum(fit)) < 2e-10),
            length(fit) == 1 + ncol(phy_gts_fixed))

  # calculate sitewise marginal likelihoods for gt=0 and gt=1
  if (length(fit) == 1L) {
    site_model_0 <- numeric(nrow(data))
    site_model_1 <- numeric(nrow(data))
    site_model_1[] <- 1
  } else {
    site_model_0 <-
      cbind(phy_gts_fixed, 0) %>% { colSums(t(.) * c(fit)) } %>%
      force_to_interval()
    site_model_1 <-
      cbind(phy_gts_fixed, 1) %>% { colSums(t(.) * c(fit)) } %>%
      force_to_interval()
  }
  site_lh_0 <- with(data, binom_likelihood(bac, dp, site_model_0, err, by_site = TRUE))
  site_lh_1 <- with(data, binom_likelihood(bac, dp, site_model_1, err, by_site = TRUE))

  # permutation test likelihoods
  perm_lhs <- map_dbl(seq_len(ncol(perm_gts)),
                      ~ sum(if_else(perm_gts[, .] == 1, site_lh_1, site_lh_0)))

  # score likelihoods for all search phylotypes
  result <-
    tibble(index = seq_len(ncol(phy_gts_search)),
           phylotype = colnames(phy_gts_search)) %>%
    mutate(lh = map_dbl(index, ~ sum(if_else(phy_gts_search[, .] == 1, site_lh_1, site_lh_0))),
           p_val = map_dbl(lh, ~ (sum(. < perm_lhs, na.rm = T) + 1) / (length(perm_lhs) + 1))) %>%
    arrange(desc(lh), p_val) %>%
    select(index, phylotype, likelihood = lh, p_val_perm = p_val)

  return(result)
}

optim_fuzzy <- function(data, phy_gts, phylo, node_dist, max_dist, mix_nodes, mix_fit) {
  # return list(optimised = logical(1), gts = interger(nrow(data)))
  neighb_nodes <-
    node_dist_range(node = mix_nodes[1],
                    phylo = phylo,
                    node_dist = node_dist,
                    max_dist = max_dist,
                    inclusive = TRUE)

  optim_sites <-
    optim_fuzzy_search(data,
                              gts_fixed = phy_gts[, mix_nodes, drop = FALSE],
                              gts_neighb = phy_gts[, neighb_nodes, drop = FALSE],
                              dist_neighb = node_dist[neighb_nodes],
                              fit = mix_fit,
                              max_sites = max_dist)

  if (length(optim_sites$n_index) == 1) {
    gt[optim_sites$site_index] %<>% { if_else(. == 0L, 1L, 0L) }
    result <- list(sites = optim_sites$site_index,
                   optimised = TRUE,
                   neighbour = neighb_nodes[optim_sites$n_index],
                   n_sites = length(optim_sites$site_index),
                   lh_delta = optim_sites$lh_delta)
  } else {
    result <- list(sites = integer(),
                   optimised = FALSE,
                   neighbour = NA_integer_,
                   n_sites = 0,
                   lh_delta = 0)
  }

  return(result)
}

#' @importFrom dplyr count mutate filter bind_cols select starts_with
#' @importFrom tidyr pivot_longer
optim_fuzzy_search <- function(data, gts_fixed, gts_neighb, dist_neighb, fit, max_sites) {
  # first column of gts_fixed corresonds to node to be optimised
  #check args
  stopifnot(is.data.frame(data),
            is.matrix(gts_fixed) && ncol(gts_fixed) > 0,
            is.matrix(gts_neighb) && ncol(gts_neighb) > 0,
            nrow(data) == nrow(gts_fixed),
            nrow(data) == nrow(gts_neighb),
            is_double(fit) && all(fit >= 0 | fit <= 1) && (abs(1 - sum(fit)) < 2e-10),
            is_scalar_integerish(max_sites))

  sites_optim <- rowSums(cbind(gts_fixed[,1], gts_neighb)) %>% { which(. > 0 & . < ncol(gts_neighb) + 1) }
  data_sub <- data[sites_optim, ]
  gts_fixed <- gts_fixed[sites_optim, , drop =F]
  gts_neighb <- gts_neighb[sites_optim, ]

  # find neighbour with highest likelihood delta after flipping up to max_dist sites
  # break ties by going with lower distance
  result <-
    tibble(site_index = sites_optim,
         gt_0 = gts_fixed[,1],
         site_lh_0 = {
           cbind(0, gts_fixed[, -1, drop = FALSE]) %>%
             { colSums(t(.) * c(fit)) } %>%
             force_to_interval() %>%
             { with(data_sub, binom_likelihood(bac, dp, ., err, by_site = TRUE)) }
         },
         site_lh_1 = {
           cbind(1, gts_fixed[, -1, drop = FALSE])  %>%
             { colSums(t(.) * c(fit)) } %>%
             force_to_interval() %>%
             { with(data_sub, binom_likelihood(bac, dp, ., err, by_site = TRUE)) }
         }) %>%
    mutate(gt_opt = if_else(site_lh_0 > site_lh_1, 0, 1)) %>%
    bind_cols(
      set_colnames(gts_neighb, str_c('N', seq_len(ncol(gts_neighb)))) %>% as_tibble()) %>%
    filter(gt_0 != gt_opt) %>%
    mutate(lh_delta = abs(site_lh_0 - site_lh_1)) %>%
    select(-gt_0, -site_lh_0, -site_lh_1) %>%
    pivot_longer(cols = starts_with('N'),
                 names_to = 'n_index',
                 names_prefix = 'N',
                 values_to = 'gt_n') %>%
    mutate(n_index = as.integer(n_index)) %>%
    filter(gt_opt == gt_n) %>%
    group_by(n_index) %>%
    arrange(desc(lh_delta)) %>%
    slice(seq_len(max_sites)) %>%
    summarise(site_index = list(site_index),
              num_alt = n(),
              lh_delta = sum(lh_delta),
              dist_n = dist_neighb[n_index[1]]) %>%
    arrange(desc(lh_delta), dist_n, n_index) %>%
    slice(1) %>%
    { `if`(nrow(.) > 0,
           with(., list(n_index = n_index, site_index = site_index[[1]], lh_delta = lh_delta)),
           list(n_index = integer(), site_index = integer(), lh_delta = 0))
    }

  return(result)
}

# Implements a binary search to optimise phylotype mixture
#' @importFrom purrr pmap map_dbl map_chr
#' @importFrom tidyr expand_grid
#' @importFrom magrittr "%T>%"
optim_phy_mix <- function(data,
                          mix_gts,
                          fit = NULL,
                          resolution = 5000L,
                          max_iter = 1000L) {

  ## TODO - only optimise over sites where nodes have different phylotypes

  stopifnot(is.data.frame(data),
            is.matrix(mix_gts),
            nrow(mix_gts) == nrow(data),
            ncol(mix_gts) > 1,
            is.null(fit) || is_double(fit) && length(fit) == ncol(mix_gts),
            is_scalar_integerish(resolution) && resolution > ncol(mix_gts),
            is_scalar_integerish(max_iter))

  if (is.null(fit)) {
    fit <- rep(1 / ncol(mix_gts), ncol(mix_gts))
  }

  sites_diff <- rowSums(mix_gts) %>% { which(. > 0 & . < ncol(mix_gts)) }
  data_sub <- data[sites_diff, ]
  t_m_gt <- t(mix_gts[sites_diff, ])

  lh_same <-
    `if`(length(sites_diff) < nrow(data),
         with(data[-sites_diff, ], binom_likelihood(bac, dp, mix_gts[-sites_diff, 1], err)),
         0)

  fit_0 <- as.integer(fit * resolution)
  resolution <- sum(fit_0)
  top_lh <-
    (colSums(t_m_gt * (fit_0 / resolution ))) %>%
    force_to_interval() %>%
    { with(data_sub, binom_likelihood(bac, dp, ., err, by_site = FALSE)) }
  top_fit <- fit_0

  n_iter <- 0L
  coeff <- 1
  search_0 <-
    expand_grid(n1 = seq_len(ncol(mix_gts)), n2 = seq_len(ncol(mix_gts))) %>%
    filter(n1 != n2)

  # memoise likelihood computations using a hasmap
  lh_hash <- hashmap::hashmap(character(), numeric())

  while(n_iter < max_iter) {
    n_iter <- n_iter + 1L

    search <-
      search_0 %>%
      mutate(fit1 = top_fit[n1],
             fit2 = top_fit[n2],
             delta = as.integer(round(coeff * fit1)),
             fit1 = fit1 - delta,
             fit2 = fit2 + delta) %>%
      filter(delta > 0)

    if (nrow(search) == 0) { break }

    top_res <-
      search %>%
      (function(x) mutate(x, fit_list = pmap(x, function(n1, n2, fit1, fit2, ...) {
        replace(top_fit, c(n1, n2), c(fit1, fit2))
      }))) %>%
      mutate(fit_hash = map_chr(fit_list, ~ digest::digest(.)),
             lh = lh_hash[[fit_hash]]) %>%
      (function(x) {
        filter(x, is.na(lh)) %>%
          mutate(lh = map_dbl(fit_list, function(fit_) {
            (colSums(t_m_gt * (fit_ / resolution ))) %>%
              force_to_interval() %>%
              { with(data_sub, binom_likelihood(bac, dp, ., err, by_site = FALSE)) }
          })) %T>% {
            with(., lh_hash[[fit_hash]] <- lh)
          }
      }) %>%
      filter(lh > top_lh) %>%
      arrange(desc(lh), n1, n2) %>%
      slice(1)

    if (nrow(top_res) > 0) {
      top_fit <- top_res$fit_list[[1]]
      top_lh <- top_res$lh[[1]]
      coeff <- 1
    } else {
      coeff <- coeff / 2
    }
  }
  return(list(prop = top_fit / resolution,
              likelihood = top_lh + lh_same,
              n_iter = n_iter))
}
