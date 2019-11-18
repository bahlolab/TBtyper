


# allele_counts should be a three dimensional array
# rows == sample, cols == variant, dim 3-1 == ref_allele count, dim 3-2 == alt_allele count
# e.g allele_counts <- array (c(ref_count, alt_count), dim c (dim(ref_count), 2))
#' @export
#' @importFrom rlang is_integer is_scalar_integerish
fit_phylotypes <- function(allele_counts,
                           phylo = get_phylo(),
                           geno = get_geno(),
                           max_phylotypes = 3L,
                           min_mix_prop = 0.005,
                           min_depth = 10L,
                           max_depth = 200L,
                           error_rate = 0.005,
                           max_p_val = 0.001,
                           n_perm = 100L,
                           seed = 0L,
                           search_span = 1L,
                           rate_perm = 0.05,
                           return_all = FALSE,
                           no_ancestors = TRUE,
                           optim_hap = TRUE
                           ) {

  # check args
  stopifnot(
    is_integer(allele_counts) && is.matrix(allele_counts) && min(allele_counts, na.rm = T) >= 0L,
    length(dim(allele_counts) == 3L) && dim(allele_counts)[3] == 2L,
    is_integer(geno) && is.matrix(geno) && max(geno, na.rm = T) == 1L && min(geno, na.rm = T) == 0L,
    is_phylo(phylo),
    setequal(rownames(geno), c(phylo$tip.label, phylo$node.label)),
    all(colnames(allele_counts) %in% colnames(geno)),
    is_scalar_integerish(max_phylotypes) & max_phylotypes > 0,
    is_scalar_integerish(min_depth) & min_depth > 0,
    is_scalar_integerish(max_depth) & max_depth > 0,
    is_scalar_proportion(error_rate),
    is_scalar_proportion(max_p_val),
    is_scalar_proportion(rate_perm))

  if (FALSE) {

    # panel <- read_rds('test_data/phylomix_panel.rds')
    # phylo <- panel$phylo
    # phylotypes <- panel$phylotypes
    # b_count <- read_rds('test_data/b_count.rds')
    # depth <- read_rds('test_data/depth.rds')
    # fun <- 'phylomatch'
    # dots <- list(search_span = Inf)
  }





  # transposed shape is more convenient for per sample testing
  b_count_t <- t(b_count)
  depth_t <- t(depth)
  phylotypes_t <- t(phylotypes)

  # ensure phylotypes_t columns are ordered by node number as colnames are dropped when adding matrix to tibble in following steps
  if (! all(colnames(phylotypes_t) == node_to_label(phylo, seq_len(Nnode2(phylo))))) {
    phylotypes_t <- phylotypes_t[, node_to_label(phylo, seq_len(Nnode2(phylo)))]
  }

  res <-
    seq_len(ncol(depth_t)) %>%
    future_map_dfr(function(i) {
      exec(phylogeneric_sample,
           phylo, b_count_t[,i], depth_t[,i], phylotypes_t,
           fun = fun, !!! dots) %>%
        mutate(sample_index = i)
    }) %>%
    select(sample_index, everything())

}


#
# # runs phylomix or phylomatch for a single sample
# phylogeneric_sample <- function(phylo, bac, dp, hts,
#                                 fun = c('phylomix', 'phylomatch'),
#                                 ) {
#   fun <- match.arg(fun)
#
#   stopifnot(
#     is_integer(bac),
#     is_integer(dp),
#     length(bac) == length(dp),
#     is_integer(hts) & is.matrix(hts),
#     nrow(hts) == length(dp),
#     is_phylo(phylo),
#     treeio::Nnode2(phylo) == ncol(hts),
#     is_scalar_proportion(err_01),
#     is_scalar_proportion(err_10),
#     is_scalar_integerish(max_phylotypes) & max_phylotypes > 1L,
#     is_scalar_integerish(min_depth) & min_depth > 0L,
#     is_scalar_integerish(max_depth) & max_depth >= min_depth,
#     is_scalar_integerish(n_perm) & n_perm >= 0L,
#     is_scalar_proportion(rate_perm) & rate_perm > 0L,
#     is_scalar_proportion(max_p_val),
#     is_scalar_proportion(min_mix_prop),
#     is_integerish(seed),
#     is_integerish(search_span) & search_span > 0L,
#     is_scalar_logical(return_all))
#
#   # TODO : better management of arguments
#   #   pass dots to argument checker
#   #   warn about ignored args
#
#   if (fun == 'phylomix') {
#     phylomix_sample(phylo, bac, dp, hts,
#                     err_01 = err_01,
#                     err_10 = err_10,
#                     max_phylotypes = max_phylotypes,
#                     min_depth = min_depth,
#                     max_depth = max_depth,
#                     max_p_val = max_p_val,
#                     min_mix_prop = min_mix_prop,
#                     search_span = search_span,
#                     optim_hap = optim_hap)
#   } else if(fun == 'phylomatch') {
#     phylomatch_sample(phylo, bac, dp, hts,
#                       min_depth = min_depth,
#                       max_depth = max_depth,
#                       err_01 = err_01,
#                       err_10 = err_10,
#                       search_span = search_span)
#   } else {
#     stop('unknown function "',  fun,'"')
#   }
#
# }

# runs phylomix for a single sample
phylomix_sample <- function(phylo, bac, dp, hts,
                            min_mix_prop = 5e-3,
                            max_p_val = 1e-3,
                            err_01 = 0.006,
                            err_10 = err_01,
                            max_phylotypes = 5L,
                            min_depth = 10L,
                            max_depth = 200L,
                            search_span = 1L,
                            optim_hap = TRUE) {

  if (FALSE) {
    panel <- read_rds('test_data/phylomix_panel.rds')
    phylo <- panel$phylo
    hts <-  t(panel$phylotypes)
    bac <- t(read_rds('test_data/b_count.rds'))[,17]
    dp <- t(read_rds('test_data/depth.rds'))[,17]
    min_mix_prop = 5e-3
    max_p_val = 1e-3
    err_01 = 0.006
    err_10 = err_01
    max_phylotypes = 5L
    min_depth = 5L
    max_depth = 200L
    search_span = 1L
  }

  data <-
    tibble(bac = bac, dp = dp, baf = bac / dp, hts = hts, ) %>%
    filter(dp >= min_depth) %>%
    mutate(dp_gt_mx = dp > max_depth,
           bac = replace(bac, dp_gt_mx, round(max_depth * baf[dp_gt_mx]) %>% as.integer()),
           dp = replace(dp, dp_gt_mx, max_depth)) %>%
    select(-dp_gt_mx) %>%
    na.omit()

  if(nrow(data) == 0) {  return(tibble(note = 'insufficient data'))  }

  hap_match <- phylomatch_sample(phylo,
                                 data = data,
                                 err_01 = err_01,
                                 err_10 = err_10,
                                 search_span = search_span)

  res <-
    filter(hap_match, (p_val <= max_p_val) | (all(p_val > max_p_val, na.rm=T))) %>%
    arrange(desc(likelihood)) %>%
    slice(1) %>%
    mutate(search = map(node, ~ filter(hap_match, node !=  .)),
           node = list(node),
           phylotype = list(phylotype),
           fit = list(1)) %>%
    select(node, phylotype, fit, likelihood, p_val, stat, df, search) %>%
    mutate(p_val = replace_na(p_val, 1))


  i <- 1L
  while((i < max_phylotypes) & (res$p_val[i] < max_p_val)) {

    # optimise previous phylotype
    if (optim_hap) {
      optim_sites <- optim_phylotype(phylo, res$node[[i]], res$fit[[i]], data, err_01, err_10)
      data$hts[optim_sites, tail(res$node[[i]], 1)] %<>% { if_else(.==0L, 1L, 0L) }
      mod <- unname(data$hts[, res$node[[i]]]) %>% { colSums(t(.) * res$fit[[i]] ) }
      res$likelihood[i] <- with(data, binom_likelihood(bac, dp, mod, err_01, err_10))
    }

    hap_mix <-
      phylomix_match_next(
        phylo,
        nodes_0 = res$node[[i]],
        lh_0 = res$likelihood[i],
        fit_0 = res$fit[[i]],
        data = data,
        err_01 = err_01,
        err_10 = err_10,
        max_p_val = max_p_val,
        search_span = search_span)

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

}

# return vector of sites to flip for better fit
# sites chosen from parents/children, whichever has the most sites
optim_phylotype <- function(phylo, nodes, fit, data, err_01, err_10) {

  node <- nodes[length(nodes)]
  nodes_0 <- setdiff(nodes, node)

  rn <- rootnode(phylo)
  adj <- c(get_children(phylo, node), parent(phylo, node))

  adj_sites <- map_df(adj, ~ tibble(nd = ., site = which(data$hts[,. ] != data$hts[, node])))

  res <-
    data[adj_sites$site, ] %>%
    mutate(site = adj_sites$site,
           ht = hts[,node],
           mod_0 = map_dbl(seq_len(n()), ~sum(c(hts[., nodes_0], 0) * fit)),
           mod_1 = map_dbl(seq_len(n()), ~sum(c(hts[., nodes_0], 1) * fit))) %>%
    mutate(mod_0 =  binom_likelihood(bac, dp, mod_0, err_01, err_10, by_site = T),
           mod_1 = binom_likelihood(bac, dp, mod_1, err_01, err_10, by_site = T)) %>%
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
                                err_01 = 0.005,
                                max_p_val = 1e-3,
                                err_10 = err_10,
                                search_span = 1L) {


  mod_0 <- data$hts[, nodes_0, drop=FALSE] %>% { colSums(t(.) * fit_0) }

  node_dist <-
    phylo %>%
    { .$edge.length %<>% replace(T, 1); .} %>%
    dist.nodes()

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
      select(baf, hts) %>%
      { .$hts <- .$hts[, node_set] ; . } %>%
      fit_mix_prop()

    mod <-
      unname(data$hts[, node_set]) %>%
      { colSums(t(.) * fit ) }

    phy_tbl$fit[[node]] <- fit
    phy_tbl$likelihood[[node]] <- with(data, binom_likelihood(bac, dp, mod, err_01, err_10))

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

phylomatch_sample <- function(phylo, bac, dp, hts,
                              data = NULL,
                              min_depth = 10L,
                              max_depth = 200L,
                              err_01 = 0.006,
                              err_10 = err_01,
                              search_span = 1L) {
  if (FALSE) {
    data <- read_rds('test_data/phylomatch_data.rds')
    bac <- data$bac
    dp <- data$dp
    hts <- data$hts
    panel <- read_rds('test_data/phylomix_panel.rds')
    phylo <- panel$phylo
    err_01 = 0.006
    err_10 = err_01
    search_span = 1L
    max_depth = 200L
    min_depth = 10L
  }

  if (is.null(data)) {
    data <-
      tibble(bac = bac, dp = dp, hts = hts, baf = bac / dp) %>%
      filter(dp >= min_depth) %>%
      mutate(dp_gt_mx = dp > max_depth,
             bac = replace(bac, dp_gt_mx, round(max_depth * baf[dp_gt_mx]) %>% as.integer()),
             dp = replace(dp, dp_gt_mx, max_depth)) %>%
      select(-dp_gt_mx, -baf) %>%
      na.omit()
  }

  data$included <- FALSE

  node_dist <-
    phylo %>%
    { .$edge.length %<>% replace(T, 1); .} %>%
    dist.nodes()

  phy_tbl <-
    as_tibble(phylo) %>%
    mutate(likelihood = 0, p_val = NA_real_, stat = NA_real_, df = NA_integer_, active = F, expanded = F)

  rn <- treeio::rootnode(phylo)
  # queue of nodes to explore
  queue <- c(rn, get_children(phylo, rn, depth = search_span))

  while(length(queue) > 0L) {
    node <- queue[1L]
    queue <- queue[-1L]
    # message('node = ', node, ' (', node_to_label(phylo, node), ')')

    children <- get_children(phylo, node, 1L)

    phy_tbl$active[node] <- T

    if (length(children) > 0L) {
      # update likelihoods
      ch_vars <- map(children, ~ which(data$hts[, node] != data$hts[, .]))
      data$included[unlist(ch_vars)] <- T

      likelihood_in <-
        map2_dbl(children, ch_vars, function(ch, vi) {
          with(data, binom_likelihood(bac[vi], dp[vi], hts[vi, ch], err_01, err_10))
        })

      likelihood_out  <-
        map_dbl(ch_vars, function(vi) {
          with(data, binom_likelihood(bac[vi], dp[vi], hts[vi, node], err_01, err_10))
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
      phy_tbl$df[node] <- with(data, sum(hts[,parent] != hts[, node]))
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
  lh_rem <- with(data, binom_likelihood(bac[rem], dp[rem], hts[rem, rn], err_01, err_10))
  phy_tbl$likelihood <- phy_tbl$likelihood + lh_rem

  # retrun table of searched nodes with likelihoods
  res <-
    filter(phy_tbl, expanded) %>%
    select(node,  phylotype = label, likelihood, p_val = p_val, stat, df) %>%
    arrange(desc(likelihood), p_val)
}

fit_mix_prop <- function(data,
                         err_01 = 0.005,
                         err_10 = 0.005,
                         correct_baf = TRUE) {

  if (FALSE) {
    method = 'baf_est'
    data <- read_rds('test_data/phylomatch_data.rds')
    data$hts <- data$hts[,c(100, 200, 300)]
  }

  stopifnot(all(c('baf', 'hts') %in% colnames(data)),
            is_double(data$baf),
            is_integer(data$hts))

  nhap <- ncol(data$hts)

  if (correct_baf) {
    data$baf %<>% { (err_10 - .) / (err_01 + err_10 - 1) }
  }

  data %>%
    filter(!rowSums(hts) %in% c(0, nhap))  %>%
    { bind_cols(select(., -hts), as_tibble(.$hts %>% set_colnames(str_c('H', seq_len(ncol(.)))))) } %>%
    mutate(site = 1:n()) %>%
    gather(starts_with('H'), key = 'ht', value = 'gt') %>%
    mutate(delta = gt - baf) %>%
    select(site, ht, delta) %>%
    spread(ht, delta) %>%
    select(-site) %>%
    as.matrix() %>%
    fit_proportions_QP()

}


fit_proportions_QP <- function(delta, check_ident = F){
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
  fit <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
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


# collapse_phylotypes <- function(gds, var_info, phylo, phylotypes, min_dist = 5L) {
#   gr <- with(var_info, GRanges(chr, IRanges(start = pos, width = 1L)))
#   seqSetFilter(gds, gr)
#
#   var_set <-
#     SeqVarTools::variantInfo(gds, expanded = TRUE) %>%
#     as_tibble() %>%
#     group_by(variant.id) %>%
#     add_count(name = 'num_allele') %>%
#     ungroup() %>%
#     { inner_join(var_info, ., c('chr', 'pos', 'ref', 'alt')) }
#
#   message('using ', nrow(var_set), ' of ', nrow(var_info), ' phylotype sites.')
#   seqSetFilter(gds, variant.id = var_set$variant.id)
#
#   phylotypes_sub <- phylotypes[, var_set$var]
#   hap_dist <- parallel_binary_dist(phylotypes_sub)
#
#   nodes_ident <-
#     which(hap_dist < min_dist & lower.tri(hap_dist), arr.ind = T) %>%
#     as_tibble() %>%
#     mutate(id = seq_len(n())) %>%
#     gather(row, col, key = 'key', value = 'node') %>%
#     select(-key) %>%
#     { left_join(., select(., node) %>% distinct() %>% mutate(depth = n_ancestor(phylo, node)), 'node')  } %>%
#     group_by(id) %>%
#     arrange(desc(depth), node) %>%
#     slice(1) %>%
#     ungroup() %>%
#     pull(node) %>%
#     unique()
#
#   phylo_sub <- `if`(length(nodes_ident) > 0L,
#                     condense_phylo(phylo, nodes_ident),
#                     phylo)
#
#   message('using ', Nnode2(phylo_sub), ' of ', Nnode2(phylo), ' phylotypes.')
#   phylotypes_sub <- phylotypes_sub[node_to_label(phylo_sub, seq_len(Nnode2(phylo_sub))),]
#
#   list(var_info = var_set, phylotypes = phylotypes_sub, phylo = phylo_sub)
#
# }
