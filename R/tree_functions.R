# test if object is phylo
is_phylo <- function(x) {

  'phylo' %in% class(x)
}

# add tbl_tree to class of dataframe
#' @importFrom rlang is_integerish is_character
as_tbl_tree <- function(df) {

  stopifnot(is.data.frame(df),
            all(c("parent", "node", "label") %in% colnames(df)),
            is_integerish(df$parent),
            is_integerish(df$node),
            is_character(df$label))

  class(df) <- union('tbl_tree', class(df))
  df
}

# convert a dataframe to phylo using tidy_tree
df2phylo <- function(df) {

  tidytree::as.phylo(as_tbl_tree(df))
}

# drop all tips from a phylo object
drop_all_tips <- function(phylo) {

  stopifnot(is_phylo(phylo))

  parent <- phylo$edge[,1]
  node <- phylo$edge[,2]
  tip_set <- node[! (node %in% parent)]
  condense_phylo(phylo, tip_set)
}

# condense phylo object removing selected nodes
#' @importFrom rlang is_integerish
#' @importFrom dplyr filter select mutate bind_rows slice
#' @importFrom magrittr "%>%"
condense_phylo <- function(phylo, nodes_rm) {

  stopifnot(is_phylo(phylo),
            is_integerish(nodes_rm))

  tbl_tree <-
    with(as_tibble(phylo), {
      for(i in nodes_rm) {
        nd = node[i]
        par <- parent[i]
        ch_i <- which(parent == nd)
        parent[ch_i] <- par
        branch.length[ch_i] <- branch.length[ch_i] + branch.length[i]
      }
      tibble(parent = parent,
             node = node,
             branch.length = branch.length,
             label = label) %>%
        slice(-nodes_rm)
    })

  new_tip <-
    filter(tbl_tree, !node %in% parent ) %>%
    pull(node)

  new_inner <- setdiff(tbl_tree$node, new_tip)

  new_order <-
    filter(tbl_tree, node %in% new_inner) %>%
    select(parent, node) %>%
    as.matrix() %>%
    t() %>% c() %>% unique()

  bind_rows(
    filter(tbl_tree, node %in% new_tip) %>%
      mutate(node = seq_along(new_tip),
             parent = length(new_tip) + match(parent, new_order)),
    filter(tbl_tree, node %in% new_inner) %>%
      mutate(parent = length(new_tip) + match(parent, new_order),
             node = length(new_tip) + match(node, new_order))) %>%
    df2phylo()
}



# attach labels to nodes in phylo describe heirachy
label_nodes_phylo <- function(phylo, sep = '.') {

  in_tree <- as_tibble(phylo)

  with(in_tree, {
    root_node <- which(parent == node)
    inner_node <- (seq_len(length(node)) > root_node)
    label[root_node] <- '0'
    label[inner_node] <- NA_character_
    i <- which((parent == root_node) & (!node == root_node) & inner_node)
    label[i] <- seq_along(i)
    for (nd in which(inner_node)) {
      which(node == nd)
      i <- which((parent == nd) & inner_node)
      label[i] <- str_c(label[match(nd, node)], seq_along(i), sep = '.')
    }
    tibble(parent = parent,
           node = node,
           branch.length = branch.length,
           label = label)
  }) %>% {
    class(.) <- c('tbl_tree', class(.));
    as.phylo(.)
  }
}



# return tibble with new labels
label_phylo <- function(phylo, node_labels,
                           sep = '|',
                           anc_sep = '-A',
                           root_label = 'root') {

  stopifnot(is_phylo(phylo),
            is_integerish(node_labels) & is_dictionaryish(node_labels),
            all(node_labels %in% phylo$edge),
            is_scalar_character(sep))


  in_tree <-
    as_tibble(phylo) %>%
    mutate(label_new = names(node_labels)[match(node, node_labels)],
           depth = n_ancestor(phylo, node),
           nchild = n_child(phylo, node),
           is_anc = (nchild == 1) & (! node %in% node_labels),
           children = map(node, function(nd) Children(phylo, nd)))

  anc <-
    filter(in_tree, is_anc) %>%
    mutate(child = map_int(node, ~ Children(phylo, .))) %>%
    with({
      level = rep(1L, length(node))
      descendent = child
      for (dp in unique(sort(setdiff(depth, max(depth)), T))) {
        ii <- which(depth == dp) %>% purrr::keep(~ child[.] %in% node)
        level[ii] <- level[match(child[ii], node)] + 1L
        descendent[ii] <- descendent[match(child[ii], node)]
      }
      tibble(label_old = label, level = level, descendent = descendent)
    })

  new_labs <-
    condense_phylo(phylo, filter(in_tree, is_anc) %>% pull(node)) %>%
    { mutate(as_tibble(.),
             depth = n_ancestor(., node),
             children = map(node, function(nd) Children(., nd)))  } %>%
    left_join(select(in_tree, label, label_new, nchild), 'label') %>%
    with({
      label_new[parent == node] <- root_label
      for (dp in seq_len(max(depth)+1L)-1L) {
        slen <- length(label_new)
        for(i in which(depth == dp)) {
          ch <- children[[i]] %>% discard(~ !is.na(label_new[.]))
          label_new[ch] <- str_c(label_new[i], sep = sep, seq_along(ch))
        }
      }
      tibble(label_old = label, label_new = label_new)
    })

  in_tree %>%
    select(parent, node, branch.length, label_old = label) %>%
    left_join(new_labs, 'label_old') %>%
    (function(df) {
      select(df, label_desc = label_new, descendent = node) %>%
        right_join(anc, 'descendent') %>%
        mutate(label_anc = str_c(label_desc, anc_sep, level)) %>%
        select(label_old, label_anc) %>%
        { left_join(df, ., 'label_old') }
    }) %>%
    mutate(label_new = if_else(is.na(label_new), label_anc, label_new)) %>%
    select(parent, node, branch.length, label = label_new) %>%
    df2phylo()
}

get_tips_by_nodelab <- function(phylo, node_lab) {
  as_tibble(phylo) %>%
    filter(label %in% node_lab) %>%
    mutate(tips = get_tip_list(phylo, node)) %>%
    with(set_names(tips, label))
}

child_lab_by_lab <- function(phylo, node_lab) {
  as_tibble(phylo) %>%
    (function(tree_tbl) {
      filter(tree_tbl, label %in% node_lab) %>%
        mutate(children = map(node, function(nd) tree_tbl$label[match(child(phylo, nd), tree_tbl$node)])) %>%
        with(set_names(children, label))
    })
}

child_dist_by_lab <- function(phylo, node_lab) {
  as_tibble(phylo) %>%
    (function(tree_tbl) {
      filter(tree_tbl, label %in% node_lab) %>%
        mutate(distance = map(node, function(nd) tree_tbl$branch.length[match(child(phylo, nd), tree_tbl$node)])) %>%
        with(set_names(distance, label))
    })
}

parent_lab_by_lab <- function(phylo, node_lab) {
  as_tibble(phylo) %>%
    (function(tree_tbl) {
      filter(tree_tbl, label %in% node_lab) %>%
        mutate(parents = map_chr(node, function(nd) tree_tbl$label[match(parent(phylo, nd), tree_tbl$node)])) %>%
        with(set_names(parents, label))
    })
}

# generic constructor for a phylo object - (ape does not provide one)
# based on ape::read.tree() and ape:::.treeBuild()
phylo_ <- function(tip_labels,
                  node_labels,
                  edge_matrix,
                  edge_lengths,
                  root_edge = NULL) {
  # rules for node numbers:
  # 1:n_tip are tips
  # ? parent is n_tip + 1 ?
  stopifnot(ncol(edge_matrix) == 2)

  phy <- list(
    edge = edge_matrix,
    edge.length = edge_lengths,
    Nnode = length(node_labels),
    node.label = node_labels,
    tip.label = tip_labels
  )
  if (!is.null(root_edge)) { phy$root.edge = root_edge }

  class(phy) <- "phylo"
  attr(phy, "order") <- "cladewise"
  phy
}

# recursive function to count ancestors
n_ancestor <- function(phylo, node) {

  parent <- phylo$edge[, 1]
  child <- phylo$edge[, 2]

  n_ancestor_rec <- function(nd) {
    i <- match(nd, child)
    if (is.na(i)) { return (0L) }
    return(n_ancestor_rec(parent[i]) + 1L)
  }
  vapply(node, n_ancestor_rec, integer(1))
}


# recursive function to number of tips bellow a node
n_offspring <- function(phylo, node) {

  parent <- phylo$edge[, 1]
  child <- phylo$edge[, 2]
  is_tip <- ! {child %in% parent}

  n_offspring_rec <- function(nd) {
    i <- which(nd == parent)
    vapply(i, function(j) {
      if (is_tip[j]) { return(1L) }
      else { return(n_offspring_rec(child[j])) }
    }, integer(1)) %>% sum()
  }
  vapply(node, n_offspring_rec, integer(1))
}

n_child <- function(phylo, node) {

  left_join(
    tibble(node = node),
    tibble(node = phylo$edge[,1]) %>%
      group_by(node) %>%
      count(),
    'node') %>%
    mutate(n = replace_na(n, 0L)) %>%
    pull(n)

}

# recursive function to number of tips bellow a node
node_n_tip <- function(phylo, node) {

  parent <- phylo$edge[, 1]
  child <- phylo$edge[, 2]
  is_tip <- ! {child %in% parent}

  node_n_tip_rec <- function(nd) {
    i <- which(nd == parent)
    vapply(i, function(j) {
      if (is_tip[j]) { return(1L) }
      else { return(m_node_n_tip_rec(child[j])) }
    }, integer(1)) %>% sum()
  }
  m_node_n_tip_rec <- memoise(node_n_tip_rec)
  vapply(node, m_node_n_tip_rec, integer(1))
}

# recursive function to get list of tips for each node
get_tip_list <- function(phylo, node, as_label=TRUE) {

  parent <- phylo$edge[, 1]
  child <- phylo$edge[, 2]
  is_tip <- ! {child %in% parent}
  label <- phylo$tip.label

  get_tip_list_rec <- function(nd) {
    i <- which(nd == parent)

    map(i, function(j) {
      if (is_tip[j]) { return(j) }
      else { return(get_tip_list_rec(child[j])) }
    }) %>% unlist()
  }
  m_get_tip_list_rec <- memoise(get_tip_list_rec)
  l <- lapply(node, m_get_tip_list_rec)

  if (as_label) {
    lapply(l, function(i) label[i])
  } else {
    lapply(l, function(i) child[i])
  }
}

label_to_node <- function(phylo, labels) {
  as_tibble(phylo) %>% with(node[match(labels, label)])
}

node_to_label <- function(phylo, nodes) {
  as_tibble(phylo) %>% with(label[match(nodes, node)])
}

get_child_list <- function(phylo, node) {

  parent <- phylo$edge[, 1]
  child <- phylo$edge[, 2]
  is_tip <- ! {child %in% parent}
  label <- phylo$tip.label

  get_child_list_rec <- function(nd) {
    i <- which(nd == parent)

    map(i, function(j) {
      if (is_tip[j]) { return(child[j]) }
      else { return(c(child[j], m_get_child_list_rec(child[j]))) }
    }) %>% unlist()
  }
  m_get_child_list_rec <- memoise(get_child_list_rec)
  map(node, m_get_child_list_rec)
}

get_children <- function(phylo, node, depth=1L) {

  stopifnot(
    is_phylo(phylo),
    is_scalar_integerish(node),
    is_scalar_integerish(depth) & depth >= 0L,
    node %in% phylo$edge)

  parent <- phylo$edge[, 1]
  children <- phylo$edge[, 2]

  get_children_rec <- function(nd, depth) {
    ch <- children[which(nd == parent)]
    if ((depth > 1L) & (length(ch) > 0L)) {
      map(ch, function(ch_) {
        get_children_rec(ch_, depth - 1L)
      }) %>% unlist() %>%
        { c(ch, .) }
    } else {
      ch
    }
  }
  if (depth > 0L) {
    get_children_rec(node, depth)
  } else {
    integer(0)
  }

}

get_child_list_by_lab <- function(phylo, labels) {
  label_to_node(phylo, labels) %>%
    get_child_list(phylo, .)%>%
    map(~ node_to_label(phylo, .)) %>%
    set_names(labels)
}

dist_to_tip <- function(phylo, node) {

  parent <- phylo$edge[, 1]
  child <- phylo$edge[, 2]
  is_tip <- ! {child %in% parent}

  dist_to_tip_rec <- function(nd) {
    i <- which(nd == parent)

    map_int(i, function(j) {
      if (is_tip[j]) { return(1L) }
      else { return(m_dist_to_tip_rec(child[j]) + 1L) }
    }) %>% min()
  }
  m_dist_to_tip_rec <- memoise(dist_to_tip_rec)
  map_int(node, m_dist_to_tip_rec)
}

get_nodes <- function(phylo) {
  phylo$edge[,1] %>% sort() %>% unique()
}

get_tips <- function(phylo) {
  phylo$edge[,2] %>% sort() %>% setdiff(get_nodes(phylo))
}

get_root <- function(phylo) {
  setdiff(unique(phylo$edge[,1]), phylo$edge[,2])
}
# recursive function to get ancestors
get_anc_list <- function(phylo, node) {

  parent <- phylo$edge[, 1]
  child <- phylo$edge[, 2]

  get_anc_rec <- function(nd) {
    i <- match(nd, child)
    if (is.na(i)) { return (integer(0)) }
    return(c(get_anc_rec(parent[i]), i))
  }
  m_get_anc_rec <- memoise(get_anc_rec)
 map(node, m_get_anc_rec)
}

get_anc_list_by_lab <- function(phylo, labels) {
  label_to_node(phylo, labels) %>%
    get_anc_list(phylo, .) %>%
    map(~ node_to_label(phylo, .))
}

ggtree_node_highlight <- function(phylo, labels) {

  as_tibble(phylo) %>%
    mutate(selected = label %in% labels,
           label = if_else(selected, label, NA_character_)) %>%
    as_tbl_tree() %>%
    as.treedata() %>% {
      ggtree(.) +
        geom_nodepoint(aes(col = label, subset = selected)) +
        geom_tippoint(aes(col = label, subset = selected)) +
        theme(legend.position = 'right')
    }
}
