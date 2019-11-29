

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
                           error_rate = 0.006,
                           max_p_val = 0.001,
                           optimise_phylotypes = FALSE,
                           optimise_fit = FALSE,
                           min_sites = 1000L
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
    is_bool(optimise_phylotypes),
    is_bool(optimise_fit))

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
      fit_sample(
        phylo = phylo_sub,
        gts = t_geno,
        sm_allele_counts = allele_counts[i, ,],
        max_phylotypes = max_phylotypes,
        min_mix_prop = min_mix_prop,
        min_depth = min_depth,
        max_depth = max_depth,
        error_rate =error_rate,
        max_p_val = max_p_val,
        optimise_phylotypes = optimise_phylotypes,
        optimise_fit = optimise_fit,
        min_sites = min_sites) %>%
        mutate(sample_id = rownames(allele_counts)[i])
    }) %>%
    select(sample_id, tidyr::everything())

  return(results)
}

#' @importFrom tidyr replace_na
#' @importFrom magrittr "%<>%"
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter if_else bind_rows tibble
fit_sample <- function(phylo,
                       gts,
                       sm_allele_counts,
                       max_phylotypes,
                       min_mix_prop,
                       min_depth,
                       max_depth,
                       error_rate,
                       max_p_val,
                       optimise_phylotypes = FALSE,
                       optimise_fit = FALSE,
                       min_sites = 1000L) {


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
  if (nrow(data) < min_sites) {  return(tibble(note = 'insufficient data'))  }

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
           fit = list(1),
           abs_diff = Inf) %>%
    select(node, phylotype, fit, likelihood, p_val, stat, df, abs_diff) %>%
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
        phylo = phylo,
        nodes_0 = res$node[[i]],
        fit_0 = res$fit[[i]],
        data = data,
        error_rate = error_rate,
        delta = min_mix_prop,
        optimise_fit = optimise_fit)

    res <-
      hap_mix %>%
      mutate(node = list(c(res$node[[i]], node)),
             phylotype = list(c(res$phylotype[[i]], phylotype))) %>%
      { bind_rows(res, .) }

    i <- i + 1L
    mod_0 <- data$gts[, res$node[[i-1]], drop=FALSE] %>% { colSums(t(.) * res$fit[[i-1]]) }
    mod_1 <- data$gts[, res$node[[i]], drop=FALSE] %>% { colSums(t(.) * res$fit[[i]]) }
    res$abs_diff[i] <- sum(abs(mod_0 - mod_1), na.rm = T)

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

#' @importFrom purrr map map2 map_dbl
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter
#' @importFrom treeio parent rootnode
#' @importFrom tidyr unnest
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

#' @importFrom purrr map map_dbl pmap_dbl
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter summarise
#' @importFrom treeio parent rootnode
#' @importFrom phangorn Ancestors Descendants
match_next <- function(phylo,
                       nodes_0,
                       fit_0,
                       data,
                       error_rate = 0.005,
                       delta = c(0.005, 0.050, 0.100),
                       exclude_ancestors = TRUE,
                       exclude_descendants = TRUE,
                       exclude_dist = 25L,
                       optimise_fit = TRUE,
                       optimise_maxiter = 500,
                       blas_threads = 2L,
                       min_rho = 0.90) {
  # note: it would be faster to store phylo as a tibble
  #   - can store meta data such as ancestors
  #   - keep snp distnace as well? (per sample)

  if (FALSE) {
    dat <- readRDS('./test/match_next.rds')
    phylo <- dat$phylo
    nodes_0 <- dat$nodes_0
    fit_0 <- dat$fit_0
    data <- dat$data
    error_rate = 0.005
    delta = c(0.005, 0.050, 0.100)
    exclude_ancestors = TRUE
    exclude_descendants = TRUE
    exclude_dist = 25L
    optimise_fit = TRUE
    optimise_maxiter = 500
    blas_threads = 2L
    min_rho = 0.90
  }

  RhpcBLASctl::blas_set_num_threads(blas_threads)
  # define nodes to exclude from search
  rn <- rootnode(phylo)
  exclude <- c(rn, nodes_0)
  if (exclude_ancestors) {
    exclude <- union(exclude, unlist(Ancestors(phylo, nodes_0, type = 'all')))
  }
  if (exclude_descendants) {
    exclude <- union(exclude, unlist(Descendants(phylo, nodes_0, type = 'all')))
  }

  node_dist <- phylo_geno_dist(phylo, magrittr::set_rownames(t(data$gts), node_to_label(phylo, seq_len(Nnode2(phylo)))))
  exclude <-
    map(nodes_0, function(n) which(node_dist[n, ] <= exclude_dist, useNames = F)) %>%
    unlist() %>%
    union(exclude)

  # return value for no candidate matches
  null_result <- tibble(node = NA_integer_, phylotype = NA_character_, likelihood = NA_real_,
                         p_val = NA_real_, stat = NA_real_, df = NA_real_, fit = list(NA_real_))

  if (length(exclude) == Nnode2(phylo)) {
    return(null_result)
  }

  # calculate marginal likelihoods under multiple hypotheses
  mod_0 <- data$gts[, nodes_0, drop=FALSE] %>% { colSums(t(.) * fit_0) }
  lh_0 <-  with(data, binom_likelihood(bac, dp, mod_0, err_01 = error_rate))

  site_lh_delta_pos <-
    vapply(delta,
           function (d) { with(data, binom_likelihood(bac, dp, pmin(mod_0 + d, 1), err_01 = error_rate, by_site = TRUE)) },
           numeric(nrow(data)))

  site_lh_delta_neg <-
    vapply(delta,
           function (d) { with(data, binom_likelihood(bac, dp, pmax(mod_0 - d, 0), err_01 = error_rate, by_site = TRUE)) },
           numeric(nrow(data)))

  result_1 <-
    as_tibble(phylo) %>%
    mutate(ancestors = map(node, ~ Ancestors(phylo, .)),
           depth = lengths(ancestors)) %>%
    mutate(data = map(node, function(node) {
      tibble(delta = delta,
             likelihood = map_dbl(seq_along(delta), function(i){
               sum(if_else(data$gts[, node] == 1, site_lh_delta_pos[,i], site_lh_delta_neg[,i]))
             }))
    })) %>%
    unnest(data) %>%
    (function(x) {
      mutate(x,
             rho = pmap_dbl(select(x, delta1 = delta, node1 = node, ancestors1 = ancestors),
                           function(delta1, node1, ancestors1) {
                             filter(x, node %in% c(node1, ancestors1), delta == delta1) %>%
                               summarise(rho = suppressWarnings(cor(depth, likelihood, method = 'spearman'))) %>%
                               pull(rho)
                           }))
    }) %>%
    filter(! node %in% exclude,
           rho >= min_rho,
           likelihood > lh_0) %>%
    select(parent, node, label, delta, likelihood) %>%
    arrange(desc(likelihood))

  if (nrow(result_1) == 0) {
    return(null_result)
  }

  # find maximal likelihood phylotype and optimise fit
  result_2 <-
    result_1[1, ] %>%
    (function(x) {
      node_set <- c(nodes_0, x$node)
      fit <- data %>% select(baf, gts) %>% { .$gts <- .$gts[, node_set] ; . } %>% fit_mix_prop()
      mod <- unname(data$gts[, node_set]) %>% { colSums(t(.) * fit ) }
      if (optimise_fit) {
        min_at <- which.min(fit)
        if (fit[min_at] > 0) {
          par0 <- fit / fit[min_at]
          bdmsk <- rep(1, length(fit)) %>% replace(min_at, 0)
          fn <- function(par) {
            par_fit <- par / sum(par)
            mod <- unname(data$gts[, node_set]) %>% { colSums(t(.) * par_fit ) }
            -with(data, binom_likelihood(bac, dp, mod, err_01 = error_rate))
          }
          opt <- suppressWarnings(
            Rcgmin::Rcgmin(par = par0, fn = fn, lower = 0.5, upper = max(par0) * 2,
                           bdmsk = bdmsk, control = list(maxit = optimise_maxiter)))
          fit <- opt$par / sum(opt$par)
        }
      }
      x %>% mutate(likelihood = with(data, binom_likelihood(bac, dp, mod, err_01 = error_rate)),
                   df = sum(mod != mod_0),
                   stat = -2 * ( lh_0 - likelihood ),
                   p_val = pchisq(stat, df, lower.tail = F),
                   fit = list(fit))
    }) %>%
    select(node, phylotype = label, likelihood, p_val = p_val, stat, df, fit)

  return(result_2)
}

# return vector of sites to flip for better fit
# sites chosen from parents/children, whichever has the most sites
#' @importFrom dplyr inner_join distinct group_by ungroup as_tibble mutate filter
#' @importFrom stringr str_remove
#' @importFrom tidyr gather
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
