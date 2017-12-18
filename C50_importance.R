
# Min-max standardization
standardize <- function(x){
(x-min(x)) * 100 / (max(x) - min(x)) 
}

# Calculate multiclass entropy (log 2)
calcEntropy <- function(x){
  pvec <- table(x)/length(x)
  log2p <- log(base = 2, ifelse(pvec == 0, 0.0000000001, pvec))
  return(-sum(pvec * log2p))
}

# Extract tree structure by recursively walking through tree
walkTree <- function(node){
  require(data.table)
  id <- node$id
  #varid <- node$split$varid
  #if(is.null(varid)) varid <- NA
  #info <- data.table(id, entropy, varid)
  
  if(is.null(node$kids)){ 
    return(NULL)
  }else{
    id_kids <- sapply(node$kids, function(x) x$id)
    info <- list(c("id" = id, id_kids))
    recursive_kids <- unlist(lapply(node$kids, walkTree), recursive = F)
    #nokids <- sapply(recursive_kids, is.null)
    # if(all(nokids)){
    #   return(as.data.frame(rbind(info)))
    # }else{
    # recursive_kids <- data.table::rbindlist(recursive_kids[!nokids], fill = TRUE)
    return(c(info, recursive_kids))
    #}
  }
}

# Calculate varimp for single tree
# Importance is the decrease in entropy for each split calculated as the difference between
# entropy before split and the weighted sum (nr of observations) of entropy in each subset
varimpC50tree <- function(tree){
  # For each of the trees get the variable id, entropy and number of observations in split
  node_info <- foreach(node = 1:length(tree), .errorhandling = "pass", .combine = rbind)%do%{
    #if(!is.null(tree$node[[node]]$kids))
    id <- node
    size <- length(tree[[node]]$fitted$`(response)`)
    entropy <- calcEntropy(tree[[node]]$fitted$`(response)`)
    variable_id <- tree[[node]]$node$split$varid
    if(is.null(variable_id)) variable_id <- NA
    
    return(data.table(id= id, entropy = entropy, size = size, variable_id = variable_id))
  }
  
  # Get the tree structure to calculate differences in entropy
  relation_info <- walkTree(tree$node)
  
  for(node in seq_along(relation_info)){
    info <- relation_info[[node]]
    parent_id <- info[1]
    children <- info[2:length(info)]
    meanDecreaseEntropy <- node_info[node_info$id == parent_id, "entropy"] - sum(node_info[node_info$id %in% children, "entropy"] * node_info[node_info$id %in% children, "size"])/node_info[node_info$id == parent_id, "size"]
    node_info[node_info$id == parent_id, "meanDecreaseEntropy"] <- meanDecreaseEntropy
  }
  
  varimpTree <- node_info[!is.na(variable_id), list("meanDecreaseEntropy" = sum(meanDecreaseEntropy, na.rm = FALSE), "splits" = .N), by = "variable_id"]
  varimpTree[, splits := splits / sum(splits)]
  
  return(varimpTree)
}

# Take C50 model, calculate varimp for each tree and average over all trees
varimpC50 <- function(model, standardize = FALSE){
  varimpList <- foreach(trial = seq_len(model[["trials"]][["Actual"]]), .inorder = TRUE) %do% {
    # Extract the trees one after another
    tree <- C50:::as.party.C5.0(model, trial = trial-1)
    return(varimpC50tree(tree))
  }
  
  varimpList <- rbindlist(varimpList)
  varimp <- varimpList[, list("meanDecreaseEntropy" = mean(meanDecreaseEntropy), "meanSplits" = mean(splits)), by = "variable_id"]
  # Match to original variable names
  varimp[, variable := model$predictors[variable_id]]  

  if(standardize == TRUE){
    set(varimp, j = "scaledMeanDecreaseEntropy", value = standardize(varimp$meanDecreaseEntropy))
  }
  
  # Bugfix https://github.com/Rdatatable/data.table/issues/939
  varimp[]
  return(varimp)
}

