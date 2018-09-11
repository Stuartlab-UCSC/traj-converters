# author: duncmc831@gmail.com
# Assume monocle_obj$State directly maps to branch assignment.
# The cell names (ids) are from sampleNames(monocle_obj@phenoData).
# Its assumed that $Pseudotime and $State and ordered the same.
#
# Graph topology is determined by reducing the
# minimumSpanningTree(monocle_object) to end nodes and branch points.

library(monocle)
library(jsonlite)

write_common_json <- function(monocle_obj, file){
  write(jsonlite::toJSON(to_common_list(monocle_obj), pretty = T), file = file)
}


write_cell_x_branch <- function(monocle_obj, file){
  write.table(to_cell_x_branch(monocle_obj),sep="\t", file=file, na = '')
}


to_common_list <- function(monocle_obj){
  graph <- make_graph(monocle_obj)
  pseudotime <- monocle_obj$Pseudotime
  branch_assignment <- monocle_obj$State
  cell_names <- colnames(monocle_obj)
  
  names(pseudotime)<- cell_names
  names(branch_assignment) <- cell_names
  
  # n<- names(branch_assignment)[1:110]
  # MST node to cell id mapping.
  id_to_mst_node <- monocle_obj@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  id_to_mst_node[,1] <- sapply(id_to_mst_node[,1], function(x){paste(c("Y_", x), collapse = "")})
  mst_node_to_ids <- list()
  for (node_name in unique(id_to_mst_node[,1])) {
    mst_node_to_ids[[node_name]] <- names(id_to_mst_node[,1][id_to_mst_node[,1] == node_name])
  }
  
  mst <- igraph::graph_from_edgelist(igraph::as_edgelist(minSpanningTree(monocle_obj)), directed =F)

  which_branch <- function(path_vertices, mst_node_to_ids, branch_assignment) {
    # returns which branch a path maps too.
    # if this path has ALL the cell ids in a branch we assume its the branch.
    total_counts <- summary(branch_assignment)
    
    # get all the cells along the path.
    cell_ids_in_path <- c()
    for (node in path_vertices){
      cell_ids_in_path <- c(cell_ids_in_path, mst_node_to_ids[[node]])  
    }
     
    count_on_branch <- summary(branch_assignment[cell_ids_in_path])
    
    branch <- which(total_counts[names(count_on_branch)] == count_on_branch, arr.ind = FALSE)
    
    indeterminant_number_of_branches <- length(branch) != 1
    if (indeterminant_number_of_branches){
      stop("The branch finding method was indeterminant.")
    }
    
    return(unname(branch)[1])
  }
  
  nodes <- igraph::V(graph)$name
  edges <- igraph::as_edgelist(graph)
  edgeIds <- c()
  nodeIds1 <- c()
  nodeIds2 <- c()
  # CM stands for cell mapping, its field in the common json format.
  cellIdCM <- c()
  edgeIdCM <- c()
  pseudotimeCM <- c()
  for (row in 1:nrow(edges)){
    node1 <- edges[row,1]
    nodeIds1 <- c(nodeIds1, node1)
    node2 <- edges[row,2]
    nodeIds2 <- c(nodeIds2, node1)
    edge_id <- paste(c(node1, node2), collapse="_")
    edgeIds <- c(edgeIds, edge_id)
    # mst nodes along the path
    path<-names(unlist(igraph::get.shortest.paths(mst,from=node1,to=node2)$vpath))
    # the branch the path represents
    branch.number <- which_branch(path, mst_node_to_ids, branch_assignment)

    n_branch_cells <- sum(branch_assignment == branch.number)
    edgeIdCM <- c(edgeIdCM, rep(edge_id, n_branch_cells))
    cellIdCM <- c(cellIdCM, names(branch_assignment[branch_assignment == branch.number]))
    pseudotimeCM <- c(pseudotimeCM, pseudotime[branch_assignment == branch.number])
  }

  output <- list(
    nodes= list(nodeId=nodes),
    egdes= list(
      edgeId=edgeIds,
      nodeId1= nodeIds1,
      nodeId2= nodeIds2
    ),
    cellMapping=list(
      cellId= cellIdCM,
      edgeId=edgeIdCM,
      psuedotime=pseudotimeCM
    )
  )
  return(output)
}


to_cell_x_branch <- function(monocle_obj){
  pseudotime <- monocle_obj$Pseudotime
  branch_assignment <- monocle_obj$State

  cell_names <- colnames(monocle_obj)
  if (length(cell_names) != length(pseudotime) || length(cell_names) != length(branch_assignment)){
    stop("Error: Are the sampleNames defined in your CellDataSet? sampleNames(monocle_obj@phenoData) is not the correct length.")
  }
  names(pseudotime)<-  cell_names
  names(branch_assignment) <- cell_names

  # TODO: implementation could be improved.
  branch_ids <- unique(branch_assignment)
  branches_list <- list()
  for (branch in branch_ids){
    cells_not_on_branch <- names(branch_assignment[branch_assignment != branch])
    cells_not_on_branch <- cells_not_on_branch[!is.na(cells_not_on_branch)]
    tmp_pseudo <- pseudotime
    tmp_pseudo[cells_not_on_branch] <- NA
    branches_list[[paste(c("branch_",branch),collapse = "")]] <- tmp_pseudo
  }
  cell_x_branch <- data.frame(branches_list, row.names = names(pseudotime))
  return(cell_x_branch)
}


make_graph <- function(monocle_obj){
  # Hack to copy the graph.
  graph <- igraph::graph_from_edgelist(igraph::as_edgelist(minSpanningTree(monocle_obj)), directed =F)
  degrees <- igraph::degree(graph)
  branches <- names(degrees[degrees>=3])
  ends <- names(degrees[degrees==1])
  middles <- names(degrees[degrees==2])
  edges_to_add <- c()
  # add the ends to the branches
  for (branch in branches){
    for (end.node in ends){
      path<-unlist(igraph::get.shortest.paths(graph,from=branch,to=end.node)$vpath)
      if (sum(!is.na(path[branches])) == 1) { # there are no other branches on that path
        edges_to_add <- c(edges_to_add, branch, end.node)
      }
    }
  }
  
  #add branch-branch edgess
  more_than_one_branch <- length(branches) > 1
  if (more_than_one_branch){
    # If the path between two branches doesn't have another branch
    # in it then connect it together.
    branch.combo <- combn(branches, 2)
    for (col in 1:ncol(branch.combo)){
      b1<- branch.combo[1,col]
      b2<- branch.combo[2,col]
      path<-unlist(igraph::get.shortest.paths(graph,from=b1,to=b2)$vpath)
      if (sum(!is.na(path[branches])) == 2) { # there are no other branches on that path
        edges_to_add <- c(edges_to_add, b1, b2)
      }
    }
  }
  
  graph <- igraph::add.edges(graph, edges_to_add)
  # Then we remove all the other nodes.
  nodes_to_remove <- middles
  graph <- igraph::delete.vertices(graph, nodes_to_remove)
  return(graph)
}