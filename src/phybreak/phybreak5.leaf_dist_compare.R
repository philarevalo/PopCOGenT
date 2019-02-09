setwd("/users/username/rumino_0_phybreak/align/")
library("ape",lib.loc ="/home/davevan/applications/R_libs/")

project_dir = "/users/username/rumino_0_phybreak/"
input_dir = paste(project_dir,"align/")
alignment_dir = paste(input_dir,"alignment_blocks/",sep="")
phy_split_dir = paste(input_dir,"phy_split/",sep="")

contig_extension = ".fa"

strain_list_file = "strain_names.txt"
output_prefix = "rumino0"

pop1 = c("NZFCFA010000011_NZFCFA010000011","ATCC29149_NZPUEL010000011","2789STDY5608852_NZCYZG010000011","RJX1120_NZNIHO010000011","RJX1121_NZNIHP010000011","CC55001C_NZKI6694141")
pop2 = c("RJX1124_NZNIHS010000011","RJX1119_NZNIHN010000011","RJX1128_NZNIHW010000011","NZFJUS010000011_NZFJUS010000011","RJX1127_NZNIHV010000011","RJX1125_NZNIHT010000011","AGR2154_NZJAGQ010000011","RJX1122_NZNIHQ010000011")

outfilename = paste(output_prefix,".core.phyml_tree_info.leaf_dists.txt", sep="")
data <- t(read.table(paste(output_prefix,".core.phyml_tree_info.txt", sep=""), header = FALSE))
trees <- as.matrix(data[6,])
mono_phy_out = paste(output_prefix,".core.phyml_tree_info.group_monophy.txt", sep="")

for(i in 1:nrow(trees))
{
  tree1 = unroot(read.tree(text=trees[i]))
  temp = cophenetic(tree1)
  pop1_mono = ""
  pop2_mono = ""
  pop1_mono = try(is.monophyletic(tree1,pop1,reroot=TRUE), silent = TRUE)
  pop2_mono = try(is.monophyletic(tree1,pop2,reroot=TRUE), silent = TRUE)
  mono = 0
  if(pop1_mono ==TRUE){
    mono = 1
  }
  if(pop2_mono ==TRUE){
    mono = 1
  }
  
  treeno = paste("##",toString(i),"##",mono, sep = "")
  write.table(treeno, file = outfilename, append = TRUE, quote = FALSE,row.names = FALSE, col.names = FALSE)
  write.table(temp, file = outfilename, append = TRUE, quote = FALSE, sep = "\t",row.names = TRUE, col.names = TRUE)
}
