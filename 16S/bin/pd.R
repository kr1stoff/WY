library(ape)

pd <- function(samp, tree, include.root = TRUE)
{
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute pd")
  }
  if (include.root) {
    if (!is.rooted(tree)) {
      stop("Rooted tree required to calculate PD with include.root=TRUE argument")
    }
    tree <- node.age(tree)
  }
  species <- colnames(samp)
  SR <- rowSums(ifelse(samp > 0, 1, 0))
  nlocations <- dim(samp)[1]
  nspecies <- dim(samp)[2]
  PDs <- rep(NA, nlocations)
  for (i in 1:nlocations) {
    present <- species[samp[i,] > 0]
    treeabsent <- tree$tip.label[which(!(tree$tip.label %in%
      present))]
    if (length(present) == 0) {
      PDs[i] <- 0
    }
    else if (length(present) == 1) {
      if (!is.rooted(tree) || !include.root) {
        warning(
          "Rooted tree and include.root=TRUE argument required to calculate PD of single-species communities. Single species community assigned PD value of NA."
        )
        PDs[i] <- NA
      }
      else {
        PDs[i] <- tree$ages[which(tree$edge[, 2] == which(tree$tip.label ==
                                                            present))]
      }
    }
    else if (length(treeabsent) == 0) {
      PDs[i] <- sum(tree$edge.length)
    }
    else {
      sub.tree <- drop.tip(tree, treeabsent)
      if (include.root) {
        if (!is.rooted(tree)) {
          stop("Rooted tree required to calculate PD with include.root=TRUE argument")
        }
        sub.tree.depth <- max(node.age(sub.tree)$ages)
        orig.tree.depth <- max(tree$ages[which(tree$edge[,
                                                 2] %in% which(tree$tip.label %in% present))])
        PDs[i] <-
          sum(sub.tree$edge.length) + (orig.tree.depth -
            sub.tree.depth)
      }
      else {
        PDs[i] <- sum(sub.tree$edge.length)
      }
    }
  }
  PDout <- data.frame(PD = PDs, SR = SR)
  rownames(PDout) <- rownames(samp)
  return(PDout)
}
