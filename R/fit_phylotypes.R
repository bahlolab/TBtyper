

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
                       intermediate_phylotypes = TRUE,
                       intermediate_max_dist = 100L) {

  # TODO - bring back optim_phylotype as slide_phylotype (could be anywhere between current node, parents and children)
  #   - give optimised a name like 1.1.2/2_0.34_1.1.2/2/1

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

  node_dist <- phylo_geno_dist(phylo, magrittr::set_rownames(t(data$gts), node_to_label(phylo, seq_len(Nnode2(phylo)))))

  # find the first phylotype
  res_1 <-
    score_phy_mix(phylo, data) %>%
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

    if (intermediate_phylotypes) {
      # move this to a separate funcition
      candidates <- node_dist_range(res_1_top$node,
                                    phylo = phylo,
                                    node_dist = node_dist,
                                    max_dist = intermediate_max_dist,
                                    inclusive = TRUE)

    }

    sample_fit <-
      res_1_top %>%
      mutate(mix_n = 1, mix_index = 1, mix_prop = 1, search = list(res_1), p_val_lrt = 0) %>%
      select(mix_n, mix_index, mix_prop, node, phylotype, mix_prop, likelihood, p_val_perm, p_val_lrt, search)

  } else {
    return(tibble(note = "no matches passing filters"))
  }

  for (i in seq_len(max_phylotypes - 1)) {
    lh_last <- filter(sample_fit, mix_n == i) %>% pull(likelihood) %>% first()
    nodes_last <- filter(sample_fit, mix_n == i) %>% pull(node)
    mix_prop_last <- filter(sample_fit, mix_n == i) %>% pull(mix_prop)

    # find next mixture component
    res_next <-
      score_phy_mix(phylo = phylo,
                    data = data,
                    nodes_fixed = nodes_last,
                    fit = c_spike_in(mix_prop_last, spike_in_p)) %>%
      arrange(desc(likelihood), p_val_perm) %>%
      filter(!node %in% exclusions(nodes_last, phylo, node_dist,
                                   exclude_parent = exclude_parent,
                                   exclude_child = exclude_child,
                                   exclude_ancestor = exclude_ancestor,
                                   exclude_descendant = exclude_descendant,
                                   exclude_distance = exclude_distance,
                                   exclude_inner = exclude_inner))

    res_next_top <-
      res_next %>%
      filter(p_val_perm < max_p_val_perm) %>% #,rho > min_rho) %>%
      slice(1)

    if (nrow(res_next_top) == 0) {
      break
    }

    mix_nodes <- c(nodes_last, res_next_top$node[1])
    optim_prop_0 <- optim_phy_mix(data, mix_nodes)
    lrt_p_val <- pchisq(-2 * (lh_last - optim_prop_0$likelihood),
                        sum(node_dist[last(mix_nodes), setdiff(mix_nodes, last(mix_nodes))]),
                        lower.tail = F)

    mix_fit <- with(optim_prop_0,
                    tibble(mix_n = length(nodes),
                           mix_index = seq_len(mix_n),
                           mix_prop = prop,
                           node = nodes,
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

#' @importFrom treeio parent rootnode
#' @importFrom phangorn Ancestors Descendants Children
#' @importFrom purrr map
exclusions <- function(nodes,
                       phylo,
                       node_dist,
                       exclude_parent = FALSE,
                       exclude_child = FALSE,
                       exclude_ancestor = FALSE,
                       exclude_descendant = FALSE,
                       exclude_distance = NULL,
                       exclude_inner = FALSE,
                       exclude_self = TRUE,
                       exclude_root = TRUE) {


  exclude <-
    c(`if`(exclude_parent, parent(phylo, nodes), integer()),
      `if`(exclude_child, unlist(Children(phylo, nodes)), integer()),
      `if`(exclude_ancestor, unlist(Ancestors(phylo, nodes)), integer()),
      `if`(exclude_descendant, unlist(Descendants(phylo, nodes)), integer()),
      `if`(exclude_root, rootnode(phylo), integer()),
      `if`(exclude_self, nodes, integer()),
      `if`(exclude_inner, inner_nodes(phylo), integer()),
      `if`(!is.null(exclude_distance),
           unlist(map(nodes, ~ which(node_dist[., ] <= exclude_distance))),
           integer())) %>%
    sort() %>%
    unique()

  return(exclude)
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
           lh = map_dbl(node, ~ sum(if_else(data$gts[, .] == 1, site_lh_1, site_lh_0))),
           p_val = map_dbl(lh, ~ (sum(. < perm_lhs, na.rm = T)) / (length(perm_lhs) + 1))) %>%
    arrange(desc(lh), p_val) %>%
    select(node, phylotype = label, likelihood = lh, p_val_perm = p_val)

  return(result)
}

optim_intermediate <- function(data, node, neighbours, max_sites) {
  # TODO - allow arbitrary states
  gt_sub <- data$gts[, c(node, neighbours)]
  sites_optim <- rowSums(gt_sub) %>% { which(. > 0 & . < ncol(gt_sub)) }
  data_sub <- select(data, bac, dp, err) %>% slice(sites_optim)
}

# Implements a binary search to optimise phylotype mixture
#' @importFrom purrr pmap map_dbl map_chr
#' @importFrom tidyr expand_grid
#' @importFrom magrittr "%T>%"
optim_phy_mix <- function(data,
                          nodes,
                          fit = NULL,
                          resolution = 5000L,
                          max_iter = 1000L) {

  ## TODO - only optimise over sites where nodes have different phylotypes

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
    expand_grid(n1 = seq_along(nodes), n2 = seq_along(nodes)) %>%
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
        # bind_rows(filter(x, !is.na(lh)),
        filter(x, is.na(lh)) %>%
          mutate(lh = map_dbl(fit_list, function(fit_) {
            (colSums(t_n_gt * (fit_ / resolution ))) %>%
              force_to_interval() %>%
              { with(data, binom_likelihood(bac, dp, ., err, by_site = FALSE)) }
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
  return(list(nodes = nodes,
              prop = top_fit / resolution,
              likelihood = top_lh,
              n_iter = n_iter))
}
