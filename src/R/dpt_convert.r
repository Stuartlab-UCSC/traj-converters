# There's a assumption here that the cell names (ids)
# can be pulled from one of the vectors in the dpt
# object instance. get.cell.names() is where this happens.
# Currently we run through each of the names in dpt object,
# and take the cell names from the first vector where names()
# is not null.

# A star topology is assumed to keep heuristics at a minimum.
# I assume the values of dpt$DPT1 are a psuedotime assignment.

write_common_json <- function(dpt_obj, file){
  write(jsonlite::toJSON(to_common_list(dpt_obj), pretty = T), file = file)
}

write_cell_x_branch <- function(dpt_obj, file){
  write.table(to_cell_x_branch(dpt_obj), file = file, sep = "\t",na = '')
}

to_cell_x_branch <- function(dpt_obj){
  #TODO: Make common module. Besides pulling what is necessary
  # out of the arg any method that has unique branch assignment
  # is going to be the same here.
  
  branch.has.value <- !is.na(dpt_obj$Branch)
  branch_assignment <- dpt_obj$Branch[branch.has.value]
  pseudotime <- dpt_obj$DPT1[branch.has.value]
  cell_names <- get.cell.names(dpt_obj)[branch.has.value]
  names(pseudotime)<-  cell_names
  names(branch_assignment) <- cell_names
  # TODO: improve implementation.
  branch_ids <- unique(branch_assignment)
  branches_list <- list()
  for (branch in branch_ids){
    cells_not_on_branch <- cell_names[branch_assignment != branch]
    tmp_pseudo <- pseudotime
    tmp_pseudo[cells_not_on_branch] <- NA
    branches_list[[paste(c("branch", branch), collapse = "_")]] <- tmp_pseudo
  }
  cell_x_branch <- data.frame(branches_list, row.names = names(pseudotime))
  return(cell_x_branch)
  
}

to_common_list <- function(dpt_obj) {
  graph <- make_graph(dpt_obj)
  branch.has.value <- !is.na(dpt_obj$Branch)
  branch_assignment <- dpt_obj$Branch[branch.has.value]
  pseudotime <- dpt_obj$DPT1[branch.has.value]
  cell_names <- get.cell.names(dpt_obj)[branch.has.value]
  names(pseudotime)<-  cell_names
  names(branch_assignment) <- cell_names
  
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
    pseudotimeCM <- c(pseudotimeCM, pseudotime[branch_assignment==branch_n])
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

make_graph <- function(dpt_obj){
  # Makes a star graph. 
  #TODO: better implementation would be use star constructor from igraph.
  branch.has.value <- !is.na(dpt_obj$Branch)
  branch_assignment <- dpt_obj$Branch[branch.has.value]
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

get.cell.names <- function(dpt_obj){
  # Look through all the vectors in dpt_obj
  # and take the first one that has names on it.
  # assuming that's the right order.
  hasNames <- names(dpt_obj)[sapply(names(dpt_obj), FUN=function(x){!is.null(names(dpt_obj[[x]]))})][1]
  return(names(dpt_obj[[hasNames]]))
}


