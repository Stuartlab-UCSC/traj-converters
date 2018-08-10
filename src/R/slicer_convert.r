

write_common_json <- function(branch_assignment, cell_ordering, cell_names, filepath){
  write(jsonlite::toJSON(to_common_list(branch_assignment, cell_ordering, cell_names), pretty = T), file = filepath)
}


write_cell_x_branch <- function(branch_assignment, cell_ordering, cell_names, filepath){
  write.table(to_cell_x_branch(branch_assignment, cell_ordering, cell_names), file=filepath, sep="\t", na='')
}


to_cell_x_branch <- function(branch_assignment, cell_ordering, cell_names){
  #TODO: Make a shared module. Keep interface matching slicer.
  pseudotime <- cell_ordering
  names(pseudotime)<-  cell_names
  names(branch_assignment) <- cell_names
  
  # TODO: improve implementation.
  branch_ids <- unique(branch_assignment)
  branches_list <- list()
  for (branch in branch_ids){
    cells_not_on_branch <- names(branch_assignment[branch_assignment != branch])
    tmp_pseudo <- pseudotime
    tmp_pseudo[cells_not_on_branch] <- NA
    branches_list[[paste(c("branch", branch), collapse = "_")]] <- tmp_pseudo
  }
  cell_x_branch <- data.frame(branches_list, row.names = names(pseudotime))
  return(cell_x_branch)
  
}


to_common_list <- function(branch_assignment, cell_ordering, cell_names){
  graph <- make_graph(branch_assignment)
  nodes <- igraph::V(graph)$name
  edges <- igraph::as_edgelist(graph)
  edgeIds <- c()
  nodeIds1 <- c()
  nodeIds2 <- c()
  # CM stands for cell mapping.
  cellIdCM <- c()
  edgeIdCM <- c()
  pseudotimeCM <- c()
  for (row in 1:(dim(edges)[1])){
    node1 <- edges[row,1]
    nodeIds1 <- c(nodeIds1, node1)
    node2 <- edges[row,2]
    nodeIds2 <- c(nodeIds2, node1)
    edge_id <- paste(c(node1,"_", node2), collapse="")
    edgeIds <- c(edgeIds, edge_id)
    branch_n <- as.numeric(tail(strsplit(node2, "_")[[1]],1))
    n_branch_cells <- sum(branch_assignment==branch_n)
    edgeIdCM <- c(edgeIdCM, rep(edge_id, n_branch_cells))
    cellIdCM <- c(cellIdCM, cell_names[branch_assignment==branch_n])
    pseudotimeCM <- c(pseudotimeCM, cell_ordering[branch_assignment==branch_n])
  }
  output <- list(
    nodes= list(nodeId=nodes),
    egdes= list(
      edgeId=edgeIds,
      nodeId1= edges[,1],
      nodeId2= edges[,2]
    ),
    cellMapping=list(
      cellId= cellIdCM,
      edgeId=edgeIdCM,
      psuedotime=pseudotimeCM
    )
  )
}


make_graph <- function(branch_assignment){
  # Makes a star graph. 
  #TODO: better implementation would be use star constructor from igraph.
  nodes <- c()
  graph <- igraph::make_empty_graph(directed=F)
  stem_name <- "stem"
  nodes <- c(nodes, stem_name)
  # Makes a star graph with node names
  graph <- igraph::add_vertices(graph,1, name=stem_name)
  for (branch in unique(branch_assignment)){
    node_name <- paste(c("end_", branch), collapse = "")
    nodes <- c(nodes, stem_name)
    graph <- igraph::add_vertices(graph, 1, name=node_name)
    graph <- igraph::add_edges(graph, c(stem_name, node_name))
  }
  return(graph)
}