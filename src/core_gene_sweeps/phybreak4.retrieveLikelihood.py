import os
import sys
import numpy

## Collect parameters
project_dir = ""
input_contig_dir = ""
contig_dir = ""
contig_extension = ""
output_prefix = ""
pop_infile_name = ""
ref_iso = ""
ref_contig = ""
focus_population = ""
len_block_threshold = 0
gap_prop_thresh = 0.0
window_size = 0
overlap = 0
MUGSY_source = ""
phyML_loc = ""
phyML_properties = ""
ape_loc = ""
percentile_threshold = 0.0
min_physplit_window_size = 0


parameter_file = open("phybreak_parameters.txt","r")
for line in parameter_file:
	line = line.strip().split(" = ")
	if len(line) > 1:
		if line[0] == "project_dir":
			project_dir = line[1].split(" #")[0]
		elif line[0] == "input_contig_dir":
			input_contig_dir = line[1].split(" #")[0]
		elif line[0] == "contig_dir":
			contig_dir = line[1].split(" #")[0]
		elif line[0] == "input_contig_extension":
			contig_extension = line[1].split(" #")[0]
		elif line[0] == "output_prefix":
			output_prefix = line[1].split(" #")[0]
		elif line[0] == "pop_infile_name":
			pop_infile_name = line[1].split(" #")[0]
		elif line[0] == "ref_iso":
			ref_iso = line[1].split(" #")[0]
		elif line[0] == "ref_contig":
			ref_contig = line[1].split(" #")[0]
		elif line[0] == "focus_population":
			focus_population = line[1].split(" #")[0]
		elif line[0] == "len_block_threshold":
			len_block_threshold = int(line[1].split(" #")[0])
		elif line[0] == "gap_prop_thresh":
			gap_prop_thresh = float(line[1].split(" #")[0])
		elif line[0] == "window_size":
			window_size = int(line[1].split(" #")[0])
		elif line[0] == "window_overlap":
			overlap = int(line[1].split(" #")[0])
		elif line[0] == "MUGSY_source":
			MUGSY_source = line[1].split(" #")[0]
		elif line[0] == "phyML_loc":
			phyML_loc = line[1].split(" #")[0]
		elif line[0] == "phyML_properties":
			phyML_properties = line[1].split(" #")[0]
		elif line[0] == "ape_loc":
			ape_loc = line[1].split(" #")[0]
		elif line[0] == "percentile_threshold":
			percentile_threshold = float(line[1].split(" #")[0])
		elif line[0] == "min_physplit_window_size":
			min_physplit_window_size = int(line[1].split(" #")[0])
parameter_file.close()

#these directories will generate if they do not already exist
input_dir = project_dir+"align/"
alignment_dir = input_dir+"alignment_blocks/"
phy_split_dir = input_dir+"phy_split/"
tree_dir = input_dir+"trees/"
msa_out_dir = input_dir+"phybreak_blocks/"

#these are output file names
strain_list_filename = "strain_names.txt"
MSA_name = output_prefix+".core.fasta"
LCB_info = output_prefix+".alignment_block_sizes.txt"
phy_prefix = output_prefix
block_loc_filename = output_prefix+".block_location.txt"
snp_loc_filename = output_prefix+".core.SNPloc.txt"
treeloc_filename = phy_prefix+".treeloc.txt"
tree_to_LCB_filename = "tree_to_LCB.txt"
lik_filename = phy_prefix+".phy_phyml_stat.txt"
Rscript_filename = "phybreak.leafdist_compare.R"
leaf_dist_file = output_prefix+".core.phyml_tree_info.leaf_dists.txt"
tree_info_file = phy_prefix+"_"+str(window_size)+".SNP_tree_summary.txt"
summary_file_name = "phybreak_result_"+focus_population+".txt"

#############   MAIN   #############
if os.path.isdir(tree_dir) == False:
	os.makedirs(tree_dir)


##Build dictionary of PhyML  likelihood values
lkfile = open(phy_split_dir+lik_filename,"r")
ML_dict = {}
subseq = phy_prefix
tree_no = ""
nonsnp = ""
lk = ""
for line in lkfile:
	line = line.strip()
	if len(line) >0:
		if line.split("[#")[0] == '. Data set ':
			tree_no = line.split("[#")[1].split("]")[0]
		elif "sites without polymorphism" in line:
			nonsnp = line.split("ism (")[1].split(").")[0]
		elif line.split(": ")[0] == '. Log likelihood of the current tree':
			temp = line.split(": ")[1]
			lk = temp[0:len(temp)-1]
			if tree_no != "":
				try:
					ML_dict[subseq][tree_no] = lk +"\t"+ nonsnp
				except:
					ML_dict[subseq] = {}
					ML_dict[subseq][tree_no] = lk +"\t"+ nonsnp
				tree_no = ""
				lk = ""
lkfile.close()

##Add location of tree in MSA to dictionary
treeloc_file = open(phy_split_dir+treeloc_filename,"r")
tree_no = 0
for line in treeloc_file:
	line = line.strip()
	tree_no +=1
	#print(tree_no)
	old = ML_dict[subseq][str(tree_no)]
	ML_dict[subseq][str(tree_no)] = line +"\t" + old

##Add tree to dictionary
tree_file = open(phy_split_dir+subseq+".phy_phyml_tree.txt","r")
tree_no = 0
for line in tree_file:
	line = line.strip()
	tree = "\t"+ line
	tree_no +=1
	ML_dict[subseq][str(tree_no)] += tree
tree_file.close()

##Write info to single file
outfile = open(input_dir+tree_info_file,"w")
subseq = phy_prefix
for tree_no in range(1,len(ML_dict[subseq])+1):
	outfile.write(ML_dict[subseq][str(tree_no)]+"\n")
	#write tree to file
	tree_num = ML_dict[subseq][str(tree_no)].split("\t")[0]
	tree_string = ML_dict[subseq][str(tree_no)].split("\t")[5]
	# tree_outfile = open(tree_dir+str(tree_num)+".nwk","w")
	# tree_outfile.write(tree_string+"\n")
	# tree_outfile.close()
outfile.close()


loc_file = open(input_dir+block_loc_filename,"r")
tree_file = open(input_dir+tree_info_file,"r")

block_start_dict = {}
block_stop_dict = {}
for line in loc_file:
	line = line.strip().split("\t")
	if line[0] != "Label":
		block_name = line[0]
		block_start_dict[block_name] = int(line[2])
		block_stop_dict[block_name] = int(line[3])
loc_file.close()

outfile = open(input_dir+tree_to_LCB_filename,"w")
outfile.write("block\tLCB\n")
for line in tree_file:
	line = line.strip().split("\t")
	tree_no = line[0]
	tree_start = int(line[1])
	tree_stop = int(line[2])
	for block_name in block_start_dict:
		if tree_start >= block_start_dict[block_name] and tree_stop <= block_stop_dict[block_name]:
			outfile.write(tree_no+"\t"+block_name+"\n")
			break
tree_file.close()
outfile.close()


### write R script for calculating pairwise distances within trees

Rscript_file = open(project_dir+Rscript_filename,"w")

Rscript_lines = 'setwd("'+input_dir+'")\n'
Rscript_lines += 'library("ape")\n'
#Make dictionary of strain to populations/groups
pop_infile = open(project_dir+pop_infile_name,"r")
focus_strain_list = []
other_strain_list = []
for line in pop_infile:
	line = line.strip().split("\t")
	strain = line[0]
	pop = line[1]
	if strain != 'Strain':
		if pop == focus_population:
			focus_strain_list.append(strain)
		else:
			other_strain_list.append(strain)
pop_infile.close()

Rscript_lines += 'pop1 = c("'
for i in range(0,len(focus_strain_list)-1):
	strain = focus_strain_list[i]
	Rscript_lines += strain +'","'
Rscript_lines += focus_strain_list[len(focus_strain_list)-1]+'")\n'

Rscript_lines += 'pop2 = c("'
for i in range(0,len(other_strain_list)-1):
	strain = other_strain_list[i]
	Rscript_lines += strain +'","'
Rscript_lines += other_strain_list[len(other_strain_list)-1]+'")\n'

Rscript_lines += 'outfilename = "'+leaf_dist_file+'"\n'
Rscript_lines += 'data <- t(read.table("'+tree_info_file+'", header = FALSE))\n'
Rscript_lines += 'trees <- as.matrix(data[6,])\n'
# Rscript_lines += 'mono_phy_out = '+ monophy_outfile+"\n"
Rscript_lines += 'for(i in 1:nrow(trees))\n{\n  tree1 = unroot(read.tree(text=trees[i]))\n  temp = cophenetic(tree1)\n  pop1_mono = ""\n  pop2_mono = ""\n  pop1_mono = try(is.monophyletic(tree1,pop1,reroot=TRUE), silent = TRUE)\n  pop2_mono = try(is.monophyletic(tree1,pop2,reroot=TRUE), silent = TRUE)\n  mono = 0\n  if(pop1_mono ==TRUE){\n    mono = 1\n  }\n  if(pop2_mono ==TRUE){\n    mono = 1\n  }\n  treeno = paste("##",toString(i),"##",mono, sep = "")\n  write.table(treeno, file = outfilename, append = TRUE, quote = FALSE,row.names = FALSE, col.names = FALSE)\n  write.table(temp, file = outfilename, append = TRUE, quote = FALSE, sep = "'+str('\t')+'",row.names = TRUE, col.names = TRUE)\n}\n'
Rscript_file.write(Rscript_lines)
Rscript_file.close()
leaf_dist_temp = open(input_dir+leaf_dist_file,"w")
leaf_dist_temp.close()
os.system("Rscript " + Rscript_filename)
