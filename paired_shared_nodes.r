# Se não funcionar...
.libPaths('/home/danimelsz/Downloads/R/x86_64-pc-linux-gnu-library/4.1')

# Remove previously loaded objects
rm(list= ls())
# Set working directory in /home/danimelsz/Desktop/Doutorado/Support_Denis/

library(ape)
library(dplyr)
library(phangorn)
library(phytools)

#### 1. FUNCTIONS ####

# Create dataframe with two columns of pairwise support values between 2 trees
pairwise_support = function(tree1, tree2){
  # Drop tips not present in both trees
  tree1 = drop.tip(tree1,setdiff(tree1$tip.label,tree2$tip.label))
  tree2 = drop.tip(tree2,setdiff(tree2$tip.label,tree1$tip.label))
  # Identify shared nodes between 2 trees
  nodes = matchNodes(tree1, tree2, method="descendants")
  nodes = as.data.frame(nodes)
  # Remove non-shared nodes
  nodes = na.omit(nodes)
  index = nodes
  # List all descendencents from each node
  d1 = Descendants(tree1, index$tr1) # descendants of tree1
  d2 = Descendants(tree2, index$tr2) # descendants of tree2
  # Note that it is required to substract the index of the internal node by the number of terminals
  nodes$tr1_updated = nodes$tr1 - tree1$Nnode - 1 # This -1 refers to the root (a node with only two edges without support value)
  nodes$tr2_updated = nodes$tr2 - tree2$Nnode - 1 
  # Create a dataframe with paired support values from 2 trees
  sn1 = nodes$tr1_updated
  sn2 = nodes$tr2_updated
  sn = data.frame(Index1 = index$tr1, Tree1 = tree1$node.label[sn1], Index2 = index$tr2, Tree2 = tree2$node.label[sn2])
  sn = sn %>% filter_all(all_vars(. != "")) # remove empty row (from the root that does not show support value)
  return(sn)
}

# Cophylo - Node Index
plot_cophylo = function(tree1, tree2, filename="cophylo_node_index.png", width=1000, height=1000, fsize=2){
  # Check if .png extension was written by the user
  if (!grepl("\\.[^.]*$", filename)) {
    filename <- paste(filename, ".png", sep = "")
  }
  # Plot
  png(filename=filename, width=width, height=height)
  #par(mfrow=c(1,1))
  #plotTree(tree1,node.numbers=T,fsize=fsize)
  #legend("topleft", legend="tree1")
  #plotTree(tree2,node.numbers=T,fsize=fsize)
  #legend("topleft", legend="tree2")
  obj = (cophylo(tree1, tree2, rotate=T, fsize=fsize,))
  plot(obj)
  nodelabels.cophylo()
  nodelabels.cophylo(which="right")
  dev.off()
}

# Cophylo - Node Support
plot_cophylo_support = function(tree1, tree2, filename="cophylo_node_support.png", width=1000, height=1000, fsize=2, cex=2){
  # Check if .png extension was written by the user
  if (!grepl("\\.[^.]*$", filename)) {
    filename <- paste(filename, ".png", sep = "")
  }
  # Plot
  png(filename=filename, width=width, height=height)
  #par(mfrow=c(1,1))
  #plotTree(tree1,node.numbers=T,fsize=fsize)
  #legend("topleft", legend="tree1")
  #plotTree(tree2,node.numbers=T,fsize=fsize)
  #legend("topleft", legend="tree2")
  obj = (cophylo(tree1, tree2, rotate=T, fsize=fsize,))
  plot(obj)
  nodelabels.cophylo(tree1$node.label,node=1:tree1$Nnode+Ntip(tree1),
                     cex=cex)
  nodelabels.cophylo(tree2$node.label,node=1:tree2$Nnode+Ntip(tree2),
                     cex=cex,which="right")
  dev.off()
}


# Plot the two trees and their strict consensus. Number in internal nodes are node indexes, not the support values
plot_2trees_1consensus = function(tree1, tree2, filename="2trees_1consensus.png", width=1000, height=5000, fsize=2){
  # Check if .png extension was written by the user
  if (!grepl("\\.[^.]*$", filename)) {
    filename <- paste(filename, ".png", sep = "")
  }
  # Plot
  png(filename=filename, width=width, height=height)
  par(mfrow=c(1,3))
  
  plotTree(tree1,node.numbers=T,fsize=fsize)
  legend("topleft", legend="tree1")
  plotTree(tree2,node.numbers=T,fsize=fsize)
  legend("topleft", legend="tree2")
  
  #plot(cophylo(tree1, tree2, rotate=F))
  
  # Consensus
  ctree = consensus(list(tree1, tree2))
  plotTree(ctree,node.numbers=T,fsize=fsize)
  legend("topleft", legend="consensus")
  dev.off()
}

#### EXAMPLE 1: SIMULATED TREES WITH 5 TIPS ####

# Read 2 fake trees
tree1 = read.newick("ver2/tree1_singleTree.nwk")
tree2 = read.newick("ver2/tree2_singleTree.nwk")

# Create dataframe with pairwise support values of shared nodes
pairwise_support(tree1, tree2)

# Plot trees and strict consensus (numbers in nodes are indexes, not support values)
plot_2trees_1consensus(tree1, tree2, filename="2trees_1consensus_example1.png")


#### EXAMPLE 2: EMPIRICAL SUBTREES WITH 9 TIPS FROM MONTESINOS (2017) ####

# Read trees
mol_hylodidae = read.newick("Support/example_trees/hylodidae_IP_mol_gb.nwk")
te_hylodidae = read.newick("Support/example_trees/hylodidae_IP_te_gb.nwk")
# Extract a clade
ex_mol = extract.clade(mol_hylodidae[[1]], 674)
ex_mol$Nnode
write.tree(ex_mol, file="ver2/clade_mol_example.nwk")
ex_te$Nnode
ex_te = extract.clade(te_hylodidae[[1]], 572)
write.tree(ex_te, file="ver2/clade_te_example.nwk")
# Plot
plot_2trees_1consensus(ex_mol, ex_te,filename="2trees_1consensus_example2.png", width=2000, height=1000, fsize=2)
plot_cophylo(ex_mol, ex_te, filename="ex2_cophylo_index.png",fsize=3, width=2000)
plot_cophylo_support(tree1=ex_mol, tree2=ex_te, filename="ex2_cophylo_support.png",fsize=3, width=2000)
# Pairwise support
data_ex2 = pairwise_support(ex_mol, ex_te)
data_ex2

