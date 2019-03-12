#!/usr/local/bin/Rscript

task <- dyncli::main()

library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(sincell)

#   ____________________________________________________________________________
#   Load data                                                               ####

expression <- as.matrix(task$expression)
params <- task$params

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# initialise sincell object
SO <- sincell::sc_InitializingSincellObject(t(expression))

# calculate distances
SO <- SO %>% sincell::sc_distanceObj(
  method = params$distance_method
)

# perform dimred, if necessary
if (params$dimred_method != "none") {
  SO <- SO %>% sincell::sc_DimensionalityReductionObj(
    method = params$dimred_method
  )
}

# cluster cells
SO <- SO %>% sincell::sc_clusterObj(
  clust.method = params$clust.method,
  mutual = params$mutual,
  max.distance = params$max.distance,
  shortest.rank.percent = params$shortest.rank.percent,
  k = params$k
)

# build graph
SO <- SO %>% sincell::sc_GraphBuilderObj(
  graph.algorithm = params$graph.algorithm,
  graph.using.cells.clustering = params$graph.using.cells.clustering,
  k = params$k_imc
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# extract cell hierarchy
cell_graph <- SO$cellstateHierarchy %>%
  igraph::as_data_frame() %>%
  rename(length = weight) %>%
  mutate(directed = FALSE)

# Leaf nodes are iteratively removed until the percentage of leaf nodes
# is below the given cutoff. Removed nodes are projected to their closest
# neighbour.
# This is to constrain the number of milestones being created.
gr <- SO$cellstateHierarchy
deg <- igraph::degree(gr)
prev_deg <- deg * 0
while (length(deg) > 10 && mean(deg <= 1) > params$pct_leaf_node_cutoff && any(deg != prev_deg)) {
  del_v <- names(which(deg == 1))
  cat("Removing ", length(del_v), " vertices with degree 1\n", sep = "")
  gr <- igraph::delete_vertices(gr, del_v)
  prev_deg <- deg[names(igraph::V(gr))]
  deg <- igraph::degree(gr)
}
to_keep <- setNames(rownames(expression) %in% names(igraph::V(gr)), rownames(expression))

#   ____________________________________________________________________________
#   Save output                                                             ####

output <- dynwrap::wrap_data(cell_ids = rownames(expression)) %>%
  dynwrap::add_cell_graph(cell_graph = cell_graph, to_keep = to_keep) %>%
  dynwrap::add_timings(timings = checkpoints)

dyncli::write_output(output, task$output)
