


# allele_counts should be a three dimensional array
# rows == sample, cols == variant, dim 3-1 == ref_allele count, dim 3-2 == alt_allele count
# e.g allele_counts <- array (c(ref_count, alt_count), dim c (dim(ref_count), 2))
#' @export
#' @importFrom rlang is_integer is_scalar_integerish is_bool
#' @importFrom treeio Nnode2
fit_phylotypes <- function(allele_counts,
                           phylo = get_phylo(),
                           geno = get_geno(),
                           max_phylotypes = 3L,
                           min_mix_prop = 0.005,
                           min_depth = 10L,
                           max_depth = 250L,
                           error_rate = 0.006,
                           max_p_val = 0.001,
                           search_span = 3L,
                           optimise_phylotypes = TRUE,
                           return_search = TRUE
                           ) {

  if (FALSE) {
    devtools::load_all()
    gds <- SeqArray::seqOpen('./test/test.gds', allow.duplicate = TRUE)
    allele_counts <- get_allele_counts_gds(gds)
    phylo = get_phylo()
    geno = get_geno()
    max_phylotypes = 3L
    min_mix_prop = 0.005
    min_depth = 10L
    max_depth = 200L
    error_rate = 0.005
    max_p_val = 0.001
    search_span = 1L
    optimise_phylotypes = TRUE
  }

  # check args
  stopifnot(
    is_integer(allele_counts) && is.array(allele_counts) && min(allele_counts, na.rm = T) >= 0L,
    length(dim(allele_counts) == 3L) && dim(allele_counts)[3] == 2L && setequal(c('Ref', 'Alt'), dimnames(allele_counts)[[3]]),
    is_integer(geno) && is.matrix(geno) && max(geno, na.rm = T) == 1L && min(geno, na.rm = T) == 0L,
    is_phylo(phylo),
    setequal(rownames(geno), c(phylo$tip.label, phylo$node.label)),
    all(colnames(allele_counts) %in% colnames(geno)),
    is_scalar_integerish(max_phylotypes) & max_phylotypes > 0,
    is_scalar_integerish(min_depth) & min_depth > 0,
    is_scalar_integerish(max_depth) & max_depth > 0,
    is_scalar_proportion(error_rate),
    is_scalar_proportion(max_p_val),
    is_bool(optimise_phylotypes))

  # subset genotypes to those in input allele counts
  geno_sub <- geno[, colnames(allele_counts)]
  message('Note: using ', ncol(geno_sub), ' of ', ncol(geno), ' genotypes')

  # subset phylotypes to those that can be distinguished over input genotypes
  phylo_sub <- collapse_phylotypes(phylo, geno_sub)
  message('Note: using ', Nnode2(phylo_sub), ' of ', Nnode2(phylo), ' phylotypes')

  # transpose genotypes
  t_geno <- t(geno_sub)[, node_to_label(phylo_sub, seq_len(Nnode2(phylo_sub)))]

  results <-
    seq_len(nrow(allele_counts)) %>%
    furrr::future_map_dfr(function(i) {
      phylomix_sample(
        phylo = phylo_sub,
        gts = t_geno,
        sm_allele_counts = allele_counts[i, ,],
        max_phylotypes = max_phylotypes,
        min_mix_prop = min_mix_prop,
        min_depth = min_depth,
        max_depth = max_depth,
        error_rate =error_rate,
        max_p_val = max_p_val,
        search_span = search_span,
        optimise_phylotypes = optimise_phylotypes) %>%
        mutate(sample_id = rownames(allele_counts)[i])
    }) %>%
    select(sample_id, tidyr::everything())

  return(results)
}

#' @importFrom tidyr replace_na
#' @importFrom magrittr "%<>%"
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter if_else bind_rows
phylomix_sample <- function(phylo,
                            gts,
                            sm_allele_counts,
                            max_phylotypes,
                            min_mix_prop,
                            min_depth,
                            max_depth,
                            error_rate,
                            max_p_val,
                            search_span,
                            optimise_phylotypes) {

  if (FALSE) {
    devtools::load_all()
    gds <- SeqArray::seqOpen('./test/test.gds', allow.duplicate = TRUE)
    allele_counts <- get_allele_counts_gds(gds)
    phylo = get_phylo()
    geno = get_geno()
    geno_sub <- geno[, colnames(allele_counts)]
    phylo_sub <- collapse_phylotypes(phylo, geno_sub)
    phylo <- phylo_sub
    t_geno <- t(geno_sub)[, node_to_label(phylo_sub, seq_len(Nnode2(phylo_sub)))]
    gts <- t_geno
    sm_allele_counts = allele_counts[1L, ,]
    max_phylotypes = 3L
    min_mix_prop = 0.005
    min_depth = 10L
    max_depth = 250L
    error_rate = 0.005
    max_p_val = 0.001
    search_span = 3L
    optimise_phylotypes = TRUE
  }

  data <-
    tibble(bac = sm_allele_counts[,'Alt'],
           dp = rowSums(sm_allele_counts),
           baf = bac / dp,
           gts = gts) %>%
    filter(dp >= min_depth) %>%
    mutate(dp_gt_mx = dp > max_depth,
           bac = replace(bac, dp_gt_mx, round(max_depth * baf[dp_gt_mx]) %>% as.integer()),
           dp = replace(dp, dp_gt_mx, max_depth)) %>%
    select(-dp_gt_mx) %>%
    na.omit()

  # arbitrary threshold on min sites
  if(nrow(data) < 10) {  return(tibble(note = 'insufficient data'))  }

  # phy_match <- phylomatch_sample(phylo,
  #                                data = data,
  #                                error_rate = error_rate,
  #                                search_span = search_span)

  phy_match <- match_first(phylo,
                           data = data,
                           error_rate = error_rate)

  res <-
    filter(phy_match, (p_val <= max_p_val) | (all(p_val > max_p_val, na.rm=T))) %>%
    arrange(desc(likelihood)) %>%
    slice(1) %>%
    mutate(search = map(node, ~ filter(phy_match, node !=  .)),
           node = list(node),
           phylotype = list(phylotype),
           fit = list(1)) %>%
    select(node, phylotype, fit, likelihood, p_val, stat, df, search) %>%
    mutate(p_val = replace_na(p_val, 1))


  i <- 1L
  while((i < max_phylotypes) & (res$p_val[i] < max_p_val)) {

    # optimise previous phylotype
    if (optimise_phylotypes) {
      optim_sites <- optim_phylotype(phylo, res$node[[i]], res$fit[[i]], data, error_rate)
      data$gts[optim_sites, tail(res$node[[i]], 1)] %<>% { if_else(.==0L, 1L, 0L) }
      mod <- unname(data$gts[, res$node[[i]]]) %>% { colSums(t(.) * res$fit[[i]] ) }
      res$likelihood[i] <- with(data, binom_likelihood(bac, dp, mod, error_rate))
    }

    hap_mix <-
      match_next(
        phylo,
        nodes_0 = res$node[[i]],
        fit_0 = res$fit[[i]],
        data = data,
        error_rate = error_rate,
        delta = min_mix_prop)

    if (min(hap_mix$fit[[1]]) < min_mix_prop) { break }

    res <-
      filter(hap_mix, (p_val <= max_p_val) | (all(p_val > max_p_val, na.rm=T))) %>%
      arrange(desc(likelihood)) %>%
      slice(1) %>%
      mutate(search = map(node, ~ filter(hap_mix, node !=  .)),
             node = list(c(res$node[[i]], node)),
             phylotype = list(c(res$phylotype[[i]], phylotype))) %>%
      { bind_rows(res, .) }

    i <- i + 1L

    # stop if we haven't found a better solution
    if (with(res, (likelihood[i] <= likelihood[i-1L]) | (is.na(p_val[i])) )) { break }
  }

  res <- arrange(res, desc(likelihood), p_val)

  return(res)
}

#' @importFrom purrr map map2_dbl map_dbl
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter
#' @importFrom treeio parent rootnode
match_first <- function(phylo,
                        data,
                        error_rate = 0.005) {


  if (FALSE) {
    error_rate = 0.005
    data <- readRDS('./test/data.rds')
    phylo <- readRDS('./test/phylo_sub.rds')
  }

  site_lh0 <- with(data, binom_likelihood(bac, dp, 0, err_01 = error_rate, by_site = TRUE))
  site_lh1 <- with(data, binom_likelihood(bac, dp, 1, err_01 = error_rate, by_site = TRUE))
  rn <- rootnode(phylo)
  rn_lh = sum(if_else(data$gts[, rn] == 1, site_lh1, site_lh0))

  result <-
    as_tibble(phylo) %>%
    mutate(data = map2(node, parent,  function(node, parent) {
      tibble(likelihood = sum(if_else(data$gts[, node] == 1, site_lh1, site_lh0))) %>%
        { `if`(node != parent,
               mutate(.,
                      df = with(data, sum(gts[, rn] != gts[, node])),
                      stat = -2 * ( rn_lh - likelihood ),
                      p_val = pchisq(stat, df, lower.tail = F)),
               .)
        }
    })) %>%
    unnest(data) %>%
    select(node,  phylotype = label, likelihood, p_val = p_val, stat, df) %>%
    arrange(desc(likelihood), p_val)

  return(result)
}

#' @importFrom purrr map map2_dbl map_dbl
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter
#' @importFrom treeio parent rootnode
#' @importFrom phangorn Ancestors Descendants
match_next <- function(phylo,
                       nodes_0,
                       fit_0,
                       data,
                       error_rate = 0.005,
                       delta = 0.005,
                       exclude_ancestors = TRUE,
                       exclude_descendants = TRUE) {



  mod_0 <- data$gts[, nodes_0, drop=FALSE] %>% { colSums(t(.) * fit_0) }
  lh_0 <-  with(data, binom_likelihood(bac, dp, mod_0, err_01 = error_rate))

  mod_0_delta_pos <- pmin(mod_0 + delta, 1)
  mod_0_delta_neg <- pmax(mod_0 - delta, 0)

  site_lh_delta_pos <- with(data, binom_likelihood(bac, dp, mod_0_delta_pos, err_01 = error_rate, by_site = TRUE))
  site_lh_delta_neg <- with(data, binom_likelihood(bac, dp, mod_0_delta_neg, err_01 = error_rate, by_site = TRUE))

  rn <- rootnode(phylo)
  exclude <- c(rn, nodes_0)
  if (exclude_ancestors) {
    exclude <- union(exclude, unlist(Ancestors(phylo, nodes_0, type = 'all')))
  }
  if (exclude_descendants) {
    exclude <- union(exclude, unlist(Descendants(phylo, nodes_0, type = 'all')))
  }

  result <-
    as_tibble(phylo) %>%
    filter(! node %in% exclude) %>%
    mutate(likelihood = map_dbl(node, function(node) {
      sum(if_else(data$gts[, node] == 1, site_lh_delta_pos, site_lh_delta_neg))
    })) %>%
    arrange(desc(likelihood)) %>%
    slice(1) %>%
    (function(x) {
      node_set <- c(nodes_0, x$node)
      fit <- data %>% select(baf, gts) %>% { .$gts <- .$gts[, node_set] ; . } %>% fit_mix_prop()
      mod <- unname(data$gts[, node_set]) %>% { colSums(t(.) * fit ) }
      x %>% mutate(likelihood = with(data, binom_likelihood(bac, dp, mod, err_01 = error_rate)),
                   df = sum(mod != mod_0),
                   stat = -2 * ( lh_0 - likelihood ),
                   p_val = pchisq(stat, df, lower.tail = F),
                   fit = list(fit))
    }) %>%
    select(node, phylotype = label, likelihood, p_val = p_val, stat, df, fit)

  return(result)
}

# return vector of sites to flip for better fit
# sites chosen from parents/children, whichever has the most sites
optim_phylotype <- function(phylo, nodes, fit, data, error_rate) {

  node <- nodes[length(nodes)]
  nodes_0 <- setdiff(nodes, node)

  rn <- treeio::rootnode(phylo)
  adj <- c(get_children(phylo, node), parent(phylo, node))

  adj_sites <- map_df(adj, ~ tibble(nd = ., site = which(data$gts[,. ] != data$gts[, node])))

  res <-
    data[adj_sites$site, ] %>%
    mutate(site = adj_sites$site,
           ht = gts[,node],
           mod_0 = map_dbl(seq_len(n()), ~sum(c(gts[., nodes_0], 0) * fit)),
           mod_1 = map_dbl(seq_len(n()), ~sum(c(gts[., nodes_0], 1) * fit))) %>%
    mutate(mod_0 =  binom_likelihood(bac, dp, mod_0, error_rate, by_site = T),
           mod_1 = binom_likelihood(bac, dp, mod_1, error_rate, by_site = T)) %>%
    select(site, ht, mod_0, mod_1) %>%
    left_join(adj_sites, 'site') %>%
    gather(mod_0, mod_1, key = 'state', value = 'lh') %>%
    mutate(state = str_remove(state, 'mod_') %>% as.integer()) %>%
    group_by(site) %>%
    arrange(desc(lh)) %>%
    slice(1) %>%
    ungroup() %>%
    filter(ht != state) %>%
    group_by(nd) %>%
    summarise(site = list(site)) %>%
    mutate(n_site = lengths(site)) %>%
    arrange(desc(n_site)) %>%
    slice(1)

  if(nrow(res) > 1) {
    unlist(res$site)
  } else{
    integer()
  }
}


phylomix_match_next <- function(phylo,
                                nodes_0,
                                lh_0,
                                fit_0,
                                data,
                                error_rate = 0.005,
                                max_p_val = 1e-3,
                                search_span = 1L) {


  mod_0 <- data$gts[, nodes_0, drop=FALSE] %>% { colSums(t(.) * fit_0) }

  node_dist <- ape::dist.nodes(phylo)

  phy_tbl <-
    as_tibble(phylo) %>%
    mutate(likelihood = 0, p_val = NA_real_, stat = NA_real_, df = NA_integer_, fit = list(numeric()))

  rn <- treeio::rootnode(phylo)
  # queue of nodes to explore
  queue <- c(rn, get_children(phylo, rn, depth = search_span))

  # remove aviod nodes from queue, add in children instead
  aq <- intersect(queue, nodes_0)
  while (length(aq > 1)) {
    ac <- map(aq, ~ get_children(phylo, .)) %>% unlist()
    queue <- setdiff(c(queue, ac), aq)
    aq <- intersect(queue, nodes_0)
  }
  searched <- integer()

  # greedy search for maximum likelihood node
  while(length(queue) > 0L) {
    node <- queue[1L]
    queue <- queue[-1L]

    node_set <- c(nodes_0, node)

    fit <-
      data %>%
      select(baf, gts) %>%
      { .$gts <- .$gts[, node_set] ; . } %>%
      fit_mix_prop()

    mod <-
      unname(data$gts[, node_set]) %>%
      { colSums(t(.) * fit ) }

    phy_tbl$fit[[node]] <- fit
    phy_tbl$likelihood[[node]] <- with(data, binom_likelihood(bac, dp, mod, error_rate))

    phy_tbl$df[node] <- sum(mod_0 != mod)
    phy_tbl$stat[node] <- -2 * ( lh_0 - phy_tbl$likelihood[node] )
    phy_tbl$p_val[node] <- pchisq(phy_tbl$stat[node], phy_tbl$df[node], lower.tail = F)

    searched <- c(searched, node)

    if (length(queue) == 0L) {
      # add additinal nodes to queue within search_span of any nodes with better likelihood than lh_0
      # and significance less than max_p_val
      queue <-
        # filter(phy_tbl, likelihood < 0, likelihood > lh_0, p_val < max_p_val) %>%
        filter(phy_tbl, likelihood < 0, likelihood > lh_0) %>%
        arrange(desc(likelihood)) %>%
        slice(1) %>%
        pull(node) %>%
        { which(node_dist[., ] <= search_span) } %>%
        setdiff(searched)
      # remove any nodes_0 nodes from queue, add their children instead
      aq <- intersect(queue, nodes_0)
      while (length(aq > 1)) {
        ac <- map(aq, ~ get_children(phylo, .)) %>% unlist()
        queue <- setdiff(c(queue, ac), aq) %>% setdiff(searched)
        aq <- intersect(queue, nodes_0)
      }
    }
  }
  # retrun table of searched nodes with likelihoods
  res <-
    phy_tbl %>%
    filter(node %in% searched) %>%
    select(node, phylotype = label, likelihood, p_val = p_val, stat, df, fit) %>%
    arrange(desc(likelihood), p_val)
}

#' @importFrom purrr map map2_dbl map_dbl
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter
#' @importFrom treeio parent rootnode
phylomatch_sample <- function(phylo,
                              data,
                              error_rate = 0.005,
                              max_p_val = 0.001,
                              search_span = 3L) {


  if (FALSE) {
    error_rate = 0.005
    max_p_val = 0.001
    search_span = 3L
    data <- readRDS('./test/data.rds')
    phylo <- readRDS('./test/phylo_sub.rds')
  }

  data$included <- FALSE

  node_dist <- ape::dist.nodes(phylo)

  phy_tbl <-
    as_tibble(phylo) %>%
    mutate(likelihood = 0, p_val = NA_real_, stat = NA_real_, df = NA_integer_, active = F, expanded = F)

  rn <- treeio::rootnode(phylo)
  # queue of nodes to explore
  # queue <- c(rn, get_children(phylo, rn, depth = search_span))
  queue <- seq_len(treeio::Nnode2(phylo))

  while(length(queue) > 0L) {
    node <- queue[1L]
    queue <- queue[-1L]
    # message('node = ', node, ' (', node_to_label(phylo, node), ')')

    children <- get_children(phylo, node, 1L)

    phy_tbl$active[node] <- T

    if (length(children) > 0L) {
      # update likelihoods
      ch_vars <- map(children, ~ which(data$gts[, node] != data$gts[, .]))
      data$included[unlist(ch_vars)] <- T

      likelihood_in <-
        map2_dbl(children, ch_vars, function(ch, vi) {
          with(data, binom_likelihood(bac[vi], dp[vi], gts[vi, ch], error_rate))
        })

      likelihood_out  <-
        map_dbl(ch_vars, function(vi) {
          with(data, binom_likelihood(bac[vi], dp[vi], gts[vi, node], error_rate))
        })

      likelihood_branch <-
        map_dbl(seq_along(children), function(ch) {
          likelihood_in[ch] + sum(likelihood_out[-ch])
        })

      phy_tbl$likelihood[children] <- { phy_tbl$likelihood[node] + likelihood_branch }
      phy_tbl$likelihood[phy_tbl$active] <- { phy_tbl$likelihood[phy_tbl$active] + sum(likelihood_out) }
    }

    phy_tbl$active[children] <- T
    phy_tbl$expanded[node] <- T

    # calculate p_val based on lr vs parent node
    if (node != rn) {
      parent <- parent(phylo, node)
      phy_tbl$df[node] <- with(data, sum(gts[,parent] != gts[, node]))
      phy_tbl$stat[node] <- -2 * ( phy_tbl$likelihood[parent] - phy_tbl$likelihood[node] )
      phy_tbl$p_val[node] <- pchisq(phy_tbl$stat[node], phy_tbl$df[node], lower.tail = F)
    }

    if (length(queue) == 0L) {
      # add additinal nodes to queue with search_span of node_ML
      node_ML <-
        filter(phy_tbl, active) %>%
        arrange(desc(likelihood)) %>%
        slice(1) %>% pull(node)
      queue <-
        which(node_dist[node_ML, ] <= search_span) %>%
        setdiff(which(phy_tbl$expanded)) %>%
        { .[order(n_ancestor(phylo, .))] }
    }
  }
  # calculate likelihood over remaining sites
  rem <- which(!data$included)
  lh_rem <- with(data, binom_likelihood(bac[rem], dp[rem], gts[rem, rn], error_rate))
  phy_tbl$likelihood <- phy_tbl$likelihood + lh_rem

  # retrun table of searched nodes with likelihoods
  result <-
    filter(phy_tbl, expanded) %>%
    select(node,  phylotype = label, likelihood, p_val = p_val, stat, df) %>%
    arrange(desc(likelihood), p_val)

  return(result)
}

