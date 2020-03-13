

# allele_counts should be a three dimensional array
# rows == sample, cols == variant, dim 3-1 == ref_allele count, dim 3-2 == alt_allele count
# e.g allele_counts <- array (c(ref_count, alt_count), dim c (dim(ref_count), 2))
#' @export
#' @importFrom rlang is_integer is_scalar_integerish is_bool
#' @importFrom treeio Nnode2
fit_phylotypes <- function(allele_counts,
                           phylo = NULL,
                           geno = NULL,
                           ref = c('h37rv', 'mrca'),
                           max_phylotypes = 3L,
                           min_mix_prop = 0.005,
                           min_depth = 50L,
                           max_depth = 250L,
                           max_p_val = 0.005,
                           optimise_phylotypes = FALSE,
                           optimise_fit = FALSE,
                           min_sites = 1000L,
                           n_perm = 999L
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
    is_scalar_proportion(max_p_val),
    is_bool(optimise_phylotypes),
    is_bool(optimise_fit),
    is_scalar_integerish(n_perm))

  # subset genotypes to those in input allele counts
  geno_sub <- geno[, colnames(allele_counts)]
  message('Note: using ', ncol(geno_sub), ' of ', ncol(geno), ' genotypes')

  # subset phylotypes to those that can be distinguished over input genotypes
  phylo_sub <- collapse_phylotypes(phylo, geno_sub)
  message('Note: using ', Nnode2(phylo_sub), ' of ', Nnode2(phylo), ' phylotypes')

  # transpose genotypes
  t_geno <- t(geno_sub)[, node_to_label(phylo_sub, seq_len(Nnode2(phylo_sub)))]
  rate <- sum(t_geno) / length(t_geno)
  perm_geno <-
    vapply(seq_len(n_perm),
           function (i) { rbinom(nrow(t_geno), 1, rate) },
           integer(nrow(t_geno)))

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
        max_p_val = max_p_val,
        min_rho = 0.95,
        spike_in_p = 0.01,
        reoptimise = TRUE,
        reoptimise_max_iter = 3L,
        min_sites = min_sites) %>%
        mutate(sample_id = rownames(allele_counts)[i])
    }) %>%
    select(sample_id, tidyr::everything())

  return(results)
}

#' @importFrom tidyr replace_na
#' @importFrom magrittr "%<>%"
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter if_else bind_rows tibble first
fit_sample <- function(phylo,
                       gts,
                       pgts,
                       sm_allele_counts,
                       max_phylotypes,
                       min_mix_prop,
                       min_depth,
                       max_depth,
                       max_p_val,
                       min_rho,
                       spike_in_p,
                       reoptimise,
                       reoptimise_max_iter,
                       min_sites) {

  # bac = b-allele count, dp = depth, gts = phylotype genotypes, pgts = permutation genotypes
  data <-
    tibble(variant = rownames(sm_allele_counts),
           bac = sm_allele_counts[,'Alt'],
           dp = rowSums(sm_allele_counts),
           baf = bac / dp,
           gts = gts,
           pgts = pgts) %>%
    filter(dp >= min_depth) %>%
    # round down samples with depth greater than max_depth
    mutate(dp_gt_mx = dp > max_depth,
           bac = replace(bac, dp_gt_mx, round(max_depth * baf[dp_gt_mx]) %>% as.integer()),
           dp = replace(dp, dp_gt_mx, max_depth)) %>%
    select(-dp_gt_mx) %>%
    # add pseudocounts and calculate 'error' rate for binomial model
    mutate(bac = bac + 1L,
           dp = dp + 2L,
           err = 1 / dp) %>%
    na.omit()

  # arbitrary threshold on min sites
  if (nrow(data) < min_sites) {  return(tibble(note = "insufficient data"))  }

  root <- rootnode(phylo)
  node_dist <- phylo_geno_dist(phylo, magrittr::set_rownames(t(data$gts), node_to_label(phylo, seq_len(Nnode2(phylo)))))

  # find the first phylotype
  res_1 <- score_phy_mix(phylo, data)
  res_1_top <-
    res_1 %>%
    filter(node != root,
           p_val < max_p_val,
           rho > min_rho) %>%
    arrange(desc(likelihood), desc(rho), p_val) %>%
    slice(1)

  # Likelihood ratio test to identify similar nodes
  res_1_lr <-
    res_1 %>%
    mutate(lr_df = node_dist[res_1_top$node, node],
           lr_stat = -2 * ( likelihood - res_1_top$likelihood),
           lr_p_val = pchisq(lr_stat, nrow(data), lower.tail = F))
  exclude_lr <- res_1_lr %>% filter(lr_p_val > 0.001) %>% pull(node)

  # store the best match in sample_fit
  if (nrow(res_1_top) > 0) {
    sample_fit <-
      res_1_top %>%
      mutate(mix_n = 1) %>%
      select(mix_n, node, phylotype, likelihood, p_val, rho) %>%
      mutate(fit = list(1), abs_diff = Inf)
  } else {
    return(tibble(note = "no matches passing filters"))
  }

  i <- 1L
  while((i < max_phylotypes) & (res$p_val[i] < max_p_val)) {

    # find next mixture component
    res_next <- score_phy_mix(phylo = phylo,
                              data = data,
                              nodes_fixed = sample_fit$node[[i]],
                              fit = c_spike_in(sample_fit$fit[[i]], spike_in_p))
    res_next_top <-
      res_next %>%
      # exclude nodes that are too closely related to the previous result
      filter(! node %in% anc_and_desc(phylo, sample_fit$node[[i]]),
             ! node %in% c(root, exclude_lr),
             p_val < max_p_val,
             rho > min_rho) %>%
      arrange(desc(likelihood), desc(rho), p_val) %>%
      slice(1)

    if (nrow(res_next_top) == 0) {
      # exit gracefully
    }

    mix_nodes <- c(sample_fit$node[[i]], res_next_top$node)
    mix_fit <- optim_phy_mix(data, mix_nodes)

    # reoptimise existing mixture components in light of new component
    if (reoptimise) {
      n_iter <- 0L
      n_alt <- 0L
      n_alt_last <- rep(-1L, i+1)
      n_alt_last[i+1] <- 0L
      while(n_iter <- reoptimise_max_iter) {
        n_iter <- n_iter + 1L
        change_made <- FALSE
        for (j in seq_along(mix_nodes)) {
          if (n_alt > n_alt_last[j]) {

          }
        }
        # exit if we couldn't find a better solution
        if (!change_made) { break }
      }
    }

    ##################################
    #  OLD STUFF
    ##################################

    res_mix <-
      match_next(
        phylo = phylo,
        nodes_last = res$node[[i]],
        fit_last = res$fit[[i]],
        data = data,
        delta = min_mix_prop,
        optimise_fit = optimise_fit,
        max_p_val = max_p_val)

    res <-
      res_mix %>%
      mutate(node = list(c(res$node[[i]], node)),
             phylotype = list(c(res$phylotype[[i]], phylotype))) %>%
      { bind_rows(res, .) }


    if (FALSE) {
      node <- first(res$node[[i]])
      phy_dist <- ape::dist.nodes(phylo)
      alts <- setdiff(which(phy_dist[node,] < optimise_radius), node)
      opt_match <- optimise_match(phylo = phylo,
                                  data = data,
                                  node = first(res$node[[i]]))
    }

    i <- i + 1L
    mod_last <- data$gts[, res$node[[i-1]], drop=FALSE] %>% { colSums(t(.) * res$fit[[i-1]]) }
    mod_1 <- data$gts[, res$node[[i]], drop=FALSE] %>% { colSums(t(.) * res$fit[[i]]) }
    res$abs_diff[i] <- sum(abs(mod_last - mod_1), na.rm = T)

    # stop if we haven't found a better solution
    if (with(res, ( (is.na(p_val[i])) ||
                    (likelihood[i] <= likelihood[i-1L]) ||
                    (min(fit[[1]]) < min_mix_prop) ))) {
      break
    }
  }

  result <- arrange(res, desc(likelihood), p_val)

  return(result)
}

#' @importFrom purrr map map_dbl pmap_dbl map2_dbl
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter summarise
#' @importFrom treeio parent rootnode
#' @importFrom phangorn Ancestors Descendants
score_phy_mix <- function(phylo, data, nodes_fixed = integer(0), fit = 1) {

  #check args
  stopifnot(is_phylo(phylo),
            is.data.frame(data),
            is_integerish(nodes_fixed),
            is_double(fit) && all(fit >= 0 | fit <= 1) && (abs(1 - sum(fit)) < 2e-10),
            length(fit) == 1 + length(nodes_fixed))

  # calculate sitewise marginal likelihoods for gt=0 and gt=1
  if (length(fit) == 1L) {
    site_model_0 <- numeric(nrow(data))
    site_model_1 <- numeric(nrow(data))
    site_model_1[] <- 1
  } else {
    site_model_0 <-
      cbind(data$gts[, nodes_fixed, drop=FALSE], 0) %>% { colSums(t(.) * c(fit)) } %>%
      force_to_interval()
    site_model_1 <-
      cbind(data$gts[, nodes_fixed, drop=FALSE], 1) %>% { colSums(t(.) * c(fit)) } %>%
      force_to_interval()
  }
  site_lh_0 <- with(data, binom_likelihood(bac, dp, site_model_0, err, by_site = TRUE))
  site_lh_1 <- with(data, binom_likelihood(bac, dp, site_model_1, err, by_site = TRUE))

  # permutation test likelihoods
  perm_lhs <- map_dbl(seq_len(ncol(data$pgts)),
                      ~ sum(if_else(data$pgts[, .] == 1, site_lh_1, site_lh_0)))

  # score likelihoods for all nodes
  result <-
    as_tibble(phylo) %>%
    mutate(ancestors = map(node, ~ Ancestors(phylo, .)),
           depth = lengths(ancestors),
           lh = map_dbl(node, ~ sum(if_else(data$gts[, .] == 1, site_lh_1, site_lh_0)))) %>%
    (function(x) {
      mutate(x, rho = map_dbl(ancestors, function(anc) {
        filter(x, node %in% anc) %>%
          with(suppressWarnings(cor(depth, lh, method = 'spearman')))
      }))
    }) %>%
    mutate(p_val = map_dbl(lh, ~ (sum(. < perm_lhs, na.rm = T) + 1) / (length(perm_lhs) + 1))) %>%
    arrange(desc(lh), desc(rho), p_val) %>%
    select(node, phylotype = label, likelihood = lh, p_val = p_val, rho)

  return(result)
}

# Implements a binary search to optimise phylotype mixture
#' @importFrom purrr pmap map_dbl
optim_phy_mix <- function(data,
                          nodes,
                          fit = NULL,
                          resolution = 5000L,
                          max_iter = 1000L) {

  stopifnot(is.data.frame(data),
            is_integer(nodes) && length(nodes) > 1,
            is.null(fit) || is_double(fit) && length(fit) == length(nodes),
            is_scalar_integerish(resolution) && resolution > length(nodes),
            is_scalar_integerish(max_iter))

  if (is.null(fit)) {
    fit <- rep(1/length(nodes), length(nodes))
  }

  t_n_gt <- t(data$gts[, nodes])

  fit_0 <- as.integer(fit * resolution)
  resolution <- sum(fit_0)
  top_lh <-
    (colSums(t_n_gt * (fit_0 / resolution ))) %>%
    force_to_interval() %>%
    { with(data, binom_likelihood(bac, dp, ., err, by_site = FALSE)) }
  top_fit <- fit_0
  fits <- list(top_fit)

  n_iter <- 0L
  coeff <- 1
  search_0 <-
    tidyr::expand_grid(n1 = seq_along(nodes), n2 = seq_along(nodes)) %>%
    filter(n1 != n2)

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
      mutate(lh = map_dbl(fit_list, function(fit_) {
        (colSums(t_n_gt * (fit_ / resolution ))) %>%
          force_to_interval() %>%
          { with(data, binom_likelihood(bac, dp, ., err, by_site = FALSE)) }
      })) %>%
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
  return(list(fit = top_fit / resolution,
              likelihood = top_lh,
              n_iter = n_iter))
}
