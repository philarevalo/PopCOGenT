##########   USER DEFINED VARIABLES   #############
import os
import sys
import numpy

project_dir = ""
contig_extension = ""
strain_list_file = ""
output_prefix = ""
pop_infile_name = ""

parameter_file = open("phybreak_parameters.txt","r")
for line in parameter_file:
	line = line.strip().split(" = ")
	if len(line) > 1:
		if line[0] == "project_dir":
			project_dir = line[1].split(" #")[0]
		elif line[0] == "contig_extension":
			contig_extension = line[1].split(" #")[0]
		elif line[0] == "strain_list_file":
			strain_list_file = line[1].split(" #")[0]
		elif line[0] == "pop_infile_name":
			pop_infile_name = line[1].split(" #")[0]
		elif line[0] == "output_prefix":
			output_prefix = line[1].split(" #")[0]
parameter_file.close()

percentile_threshold = 0.05
min_physplit_window_size = 10 #minimum number of SNPs in a row that satisfy the monophyly and percentile threshold to warrant creating a new range to output

phyML_loc = "PhyML "
phyML_properties = "-q -m JC69 -f e -c 2 -a 0.022" #-m JC69 -f e -c $Ncat -a $alpha 
slurm_prefix = "#!/bin/bash\n#SBATCH -N 1#SBATCH -n 1\n#SBATCH -p sched_mit_chisholm,sched_mit_hill,newnodes\n#SBATCH --mem=2000\n#SBATCH --time=48:00:00\ncd "+project_dir+"\n"


#############   FUNCTIONS   #############

def tree_dist_sum(tree_in):
	tree_in = tree_in.replace(",(","$").replace(",","$").replace(")","$").split("$")
	sum_out = 0.0
	for i in range(0,len(tree_in)):
		if ":" in tree_in[i]:
			num = float(tree_in[i].split(":")[1])
			sum_out += num
	return sum_out

def msa_subset(seq_dict,strt,stp):
	seq_out = {}
	for iso in seq_dict:
		seq = seq_dict[iso][strt:stp+1]
		seq_out[iso] = seq
	return seq_out

#fills in space for phylip format sequence headers
def space_fill(text,length):
	num_spaces = length-len(text)
	spaces = ""
	for i in range(0,num_spaces):
		spaces += " "
	text_out = text+spaces
	return text_out

#writes MSA info in PHYLIP format
def fasta_2_phylip(seqs_dict,window_size):
	outseq = " "+str(len(seqs_dict))+" "+str(window_size)+"\n"
	head_dict = {}
	a=-1
	for header in seqs_dict:
		a += 1
		head_dict[a] = header
	for i in range(0,len(head_dict)):
		header = head_dict[i]
		outseq += space_fill(header,10)
		outseq += " "+ seqs_dict[header] + "\n"
	outseq += "\n"
	return outseq

def dict_to_fasta(seq_dict_in):
	out_string = ''
	for iso in seq_dict_in:
		out_string += ">"+iso+"\n"+seq_dict_in[iso]+"\n"
	return out_string

##################### MAIN #####################
input_dir = project_dir+"align/"
alignment_dir = input_dir+"alignment_blocks/"
phy_split_dir = input_dir+"phy_split/"
msa_out_dir = input_dir+"phybreak_blocks/"

MSA_name = output_prefix+".core.fasta"
phy_prefix = output_prefix+".core"
trees_dist  =  output_prefix+".core.phyml_tree_info.leaf_dists.txt"
summary_file_name = output_prefix+".core.phyml_tree_info.leaf_dists.group_dist.txt"
mono_phy_name = output_prefix+".core.phyml_tree_info.group_monophy.txt"

##Make dictionary of strain to populations/groups
pop_infile = open(project_dir+pop_infile_name,"r")
pop_dict = {}
pop_list = []
for line in pop_infile:
	line = line.strip().split("\t")
	strain = line[0]
	pop = line[1]
	pop_dict[strain] = pop
	pop_list.append(pop)
pop_infile.close()
pop_list = list(set(pop_list))

##collate the distances in the leaf_dist output file

dist_dict = {}
max_dist_dict = {}
mono_phy_dict = {}
tree_no = 0
last_tree_no = 0
a = 0
strain_list = []
used = {}

all_pairwise_sum = {}
infile = open(input_dir+trees_dist,"r")
for line in infile:
	line = line.strip()
	if line[0] == "#":
		tree_no = line.split("##")[1]
		monophy = line.split("##")[2]
		last_tree_no = int(tree_no)
		used = {}
		dist_dict[tree_no] = {}
		mono_phy_dict[tree_no] = monophy
		max_dist_dict[tree_no] = 0.0
		all_pairwise_sum[tree_no] = 0.0
		a = 1 #next line is horizontal list of strain names
	elif a == 1:
		strain_list = line.split("\t")
		a = 0
	else:
		iso1 = line.split("\t")[0]
		dists = line.split("\t")
		for i in range(1,len(dists)):
			iso2 = strain_list[i-1]
			if iso1 != iso2:
				pop1 = pop_dict[iso1]
				pop2 = pop_dict[iso2]
				pairF = iso1+"\t"+iso2
				pairR = iso2+"\t"+iso1
				pop_pairF = pop1+"\t"+pop2
				try:
					used[pairF]
				except:
					try:
						used[pairR]
					except:
						dist = float(dists[i])#/branch_sum_dict[tree_no]
						all_pairwise_sum[tree_no] += dist
						used[pairF] = ""
						used[pairR] = ""
						try:
							dist_dict[tree_no][pop_pairF].append(dist)
						except:
							dist_dict[tree_no][pop_pairF] = []
							dist_dict[tree_no][pop_pairF].append(dist)
						
						if dist > max_dist_dict[tree_no]:
							max_dist_dict[tree_no] = dist
infile.close()
del used
print("done iterating over leaf_dist file")

##collect distances to calculate p-value
dict_collect = {}
for tree_no in dist_dict:
	for i in range(0,len(pop_list)):
		pop1 = pop_list[i]
		pop_pair = pop1+"\t"+pop1
		dist_list = dist_dict[tree_no][pop_pair]
		# rel_dist_list = []
		# for num in dist_list:
		# 	rel_dist_list.append(num)
		# avg = str(numpy.average(dist_list))
		dist_sum = str(numpy.sum(dist_list)/all_pairwise_sum[tree_no])
		try:
			dict_collect[pop_pair].append(dist_sum)
		except:
			dict_collect[pop_pair] = []
			dict_collect[pop_pair].append(dist_sum)
print("done collecting distances")

##Search the distances for max, min, and average distance within populations for each locus
percentile_dict = {}
physplit_trees = []
max_dist_counted = {}
outfile = open(input_dir+summary_file_name,"w")
outfile.write("tree_no\tmonophy")
for i in range(0,len(pop_list)):
	pop = pop_list[i]
	max_dist_counted[pop] = 0.0
	outfile.write("\t"+pop+"\t")
outfile.write("\n")
dist_list = []
for j in range(1,last_tree_no+1): # in dist_dict:
	tree_no = str(j)
	outfile.write(tree_no+"\t"+mono_phy_dict[tree_no])
	min_percentile = 999.9
	for i in range(0,len(pop_list)):
		pop1 = pop_list[i]
		pop_pair = pop1+"\t"+pop1
		dist_list = dist_dict[tree_no][pop_pair]
		pop_dist_list = dict_collect[pop_pair]
		dist_sum = str(numpy.sum(dist_list)/all_pairwise_sum[tree_no])
		percentile = float(sum(x <= dist_sum for x in pop_dist_list))/float(len(pop_dist_list))
		outfile.write("\t"+dist_sum+"\t"+str(percentile))
		if percentile < min_percentile:
			min_percentile = percentile
			max_dist = (dist_sum,pop1)
	percentile_dict[tree_no] = min_percentile
	if min_percentile <= percentile_threshold and mono_phy_dict[tree_no] == "1":
		physplit_trees.append(tree_no)
		outfile.write("\t1")
		if max_dist[0] > max_dist_counted[max_dist[1]]:
			max_dist_counted[max_dist[1]] = max_dist[0]
	else:
		outfile.write("\t0")
	outfile.write("\n")
outfile.close()