rm(list = ls(all.names = TRUE))

library(phytools)
library(ape)

# Export this function just to get the editted midpoint rooting function working.
getAncestors<-function(tree,node,type=c("all","parent")){
  if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
  type<-type[1]
  if(type=="all"){
    aa<-vector()
    rt<-Ntip(tree)+1
    currnode<-node
    while(currnode!=rt){
      currnode<-getAncestors(tree,currnode,"parent")
      aa<-c(aa,currnode)
    }
    return(aa)
  } else if(type=="parent"){
    aa<-tree$edge[which(tree$edge[,2]==node),1]
    return(aa)
  } else stop("do not recognize type")
}

# *Slightly* modified function that is different only by one line.
midpoint.root_EDIT<-function(tree){
  D<-cophenetic(tree)
  dd<-max(D)
  ii<-which(D==dd)[1]
  ii<-c(ceiling(ii/nrow(D)),ii%%nrow(D))
  if(ii[2]==0) ii[2]<-nrow(D)
  spp<-rownames(D)[ii]
  nn<-which(tree$tip.label==spp[2])
  tree<-reroot(tree,nn,tree$edge.length[which(tree$edge[,2]==nn)])
  aa<-getAncestors(tree,which(tree$tip.label==spp[1]))
  D<-c(0,dist.nodes(tree)[which(tree$tip.label==spp[1]),aa])
  names(D)[1]<-which(tree$tip.label==spp[1])
  i<-0
  while(D[i+1]<(dd/2)) i<-i+1
  if(i == 0) { i <- 1} # This is the only added line, just to avoid an error with some really unusual gene trees.
  tree<-reroot(tree,as.numeric(names(D)[i]),D[i+1]-dd/2)
  tree
}

# Change to this function to avoid scientific notation in output trees.
# Based on: https://www.mail-archive.com/r-sig-phylo@r-project.org/msg03045.html
# .write.tree2 <- function(phy, digits = 10, tree.prefix = "", check_tips)
# {
#   brl <- !is.null(phy$edge.length)
#   nodelab <- !is.null(phy$node.label)
#   if (check_tips) phy$tip.label <- checkLabel(phy$tip.label)
#   if (nodelab) phy$node.label <- checkLabel(phy$node.label)
#   #### NOTE -- this is the line I changed:
#   f.d <- paste(":%.", digits, "f", sep = "")
#   ### f.d <- paste0(":%.", digits, "g")
#   n <- length(phy$tip.label)
#   
#   ## terminal branches:
#   terms <- match(seq_len(n), phy$edge[, 2])
#   TERMS <- phy$tip.label
#   if (brl) TERMS <- paste0(TERMS, sprintf(f.d, phy$edge.length[terms]))
#   
#   ## internal branches, including root edge:
#   INTS <- rep(")", phy$Nnode)
#   if (nodelab) INTS <- paste0(INTS, phy$node.label)
#   if (brl) {
#     tmp <- phy$edge.length[-terms][order(phy$edge[-terms, 2])]
#     tmp <- c("", sprintf(f.d, tmp))
#     if (!is.null(phy$root.edge)) tmp[1L] <- sprintf(f.d, phy$root.edge)
#     INTS <- paste0(INTS, tmp)
#   }
#   ###    INTS[1] <- paste0(INTS[1], ";")
#   
#   ###    ## borrowed from phangorn:
#   ###    parent <- phy$edge[, 1]
#   ###    children <- phy$edge[, 2]
#   ###    kids <- vector("list", n + phy$Nnode)
#   ###    for (i in 1:length(parent))
#   ###        kids[[parent[i]]] <- c(kids[[parent[i]]], children[i])
#   ###    Nkids <- lengths(kids, FALSE)
#   ###    root <- parent[! parent %in% children][1]
#   
#   ## find the root node:
#   tmp.nodes <- unique.default(phy$edge[, 1L])
#   tmp.m <- match(tmp.nodes, phy$edge[, 2L])
#   root <- tmp.nodes[is.na(tmp.m)]
#   if (length(root) > 1)
#     stop("seems there is more than one root node")
#   storage.mode(root) <- "integer"
#   
#   ###    o <- postorder(phy)
#   o <- reorderRcpp(phy$edge, n, root, 2L)
#   ANC <- phy$edge[o, 1L]
#   DESC <- phy$edge[o, 2L]
#   NEWICK <- character(n + phy$Nnode)
#   NEWICK[1:n] <- TERMS
#   ###    root <- n + 1L
#   from <- to <- 1L
#   repeat {
#     thenode <- ANC[from]
#     if (thenode == root) {
#       to <- length(ANC)
#     } else {
#       while (ANC[to + 1L] == thenode) to <- to + 1L
#     }
#     tmp <- paste(NEWICK[DESC[from:to]], collapse = ",")
#     tmp <- paste0("(", tmp, INTS[thenode - n])
#     NEWICK[thenode] <- tmp
#     if (thenode == root) break
#     from <- to + 1L
#   }
#   paste0(NEWICK[root], ";")
# }
# assignInNamespace(".write.tree2", .write.tree2, "ape")
###assignInNamespace("midpoint.root", midpoint.root_EDIT, "phytools")


args = commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Exactly four input arguments needed: input and output tree filenames, TRUE/FALSE to indicate whether alphanumeric characters should be removed from tip labels, and TRUE/FALSE to indicate whether edge lengths should be removed..",
       call. = FALSE)
}

input_tree <- ape::read.tree(args[1])

if (! ape::is.binary(input_tree)) {
  input_tree <- ape::multi2di(input_tree) 
}

if (! ape::is.rooted(input_tree)) {
  input_tree <- midpoint.root_EDIT(input_tree)
}

if (as.logical(args[3])) {
  input_tree$tip.label <- gsub("_", "", input_tree$tip.label)
  input_tree$tip.label <- gsub("-", "", input_tree$tip.label)
  input_tree$tip.label <- gsub("\\.", "", input_tree$tip.label)
  input_tree$tip.label <- gsub(",", "", input_tree$tip.label)
  input_tree$tip.label <- make.unique(input_tree$tip.label)
}

input_tree$node.label <- NULL

if (as.logical(args[4])) {
  input_tree$edge.length <- NULL
}

ape::write.tree(phy = input_tree, file = args[2])
