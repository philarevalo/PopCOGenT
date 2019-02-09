args <- commandArgs(trailingOnly = TRUE)
setwd("/nobackup1b/users/davevan/pop_genomes/ruminococcus/rumino_0_phybreak/align/")
library("ape",lib.loc ="/home/davevan/applications/R_libs/")
output_prefix = "rumino0"

rand_itt_num = args[1]
green = strsplit(args[2],split=",")[[1]]
red = strsplit(args[3],split=",")[[1]]

data <- t(read.table(paste(output_prefix,".core.phyml_tree_info.txt", sep=""), header = FALSE))
trees <- as.matrix(data[6,])
mono_phy_out = paste(output_prefix,".rand_itt_",rand_itt_num,".group_monophy.txt", sep="")

for(i in 1:nrow(trees))
{
  tree1 = unroot(read.tree(text=trees[i]))
  green_mono = try(is.monophyletic(tree1,green,reroot=TRUE), silent = TRUE)
  red_mono = try(is.monophyletic(tree1,red,reroot=TRUE), silent = TRUE)
  mono = '0'
  treeno = paste("##",toString(i), sep = "")
  if(green_mono ==TRUE){
    mono = '1'
  }
  if(red_mono ==TRUE){
    mono = '1'
  }
  write.table(paste(toString(i),mono,sep="\t"),file=mono_phy_out, append = TRUE, quote = FALSE,row.names = FALSE, col.names = FALSE)
}
