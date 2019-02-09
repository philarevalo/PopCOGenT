##########   USER DEFINED VARIABLES   #############
import os
import numpy
from random import randint

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
		elif line[0] == "output_prefix":
			output_prefix = line[1].split(" #")[0]
		elif line[0] == "pop_infile_name":
			pop_infile_name = line[1].split(" #")[0]
		# elif line[0] == "contig_dir":
		# 	contig_dir = line[1].split(" #")[0]
		# elif line[0] == "ref_iso":
		# 	ref_iso = line[1].split(" #")[0]
		# elif line[0] == "ref_contig":
		# 	ref_contig = line[1].split(" #")[0]
		# elif line[0] == "len_block_threshold":
		# 	len_block_threshold = int(line[1].split(" #")[0])
		# elif line[0] == "gap_prop_thresh":
		# 	gap_prop_thresh = float(line[1].split(" #")[0])
		# elif line[0] == "window_size":
		# 	window_size = int(line[1].split(" #")[0])
parameter_file.close()

input_dir = project_dir+"align/"
alignment_dir = input_dir+"alignment_blocks/"
phy_split_dir = input_dir+"phy_split/"

MSA_name = output_prefix+".core.fasta"
phy_prefix = output_prefix+".core"

#region = 'pBlocks'
region = 'core'
trees_dist  =  output_prefix+".core.phyml_tree_info.leaf_dists.txt"
outfile_name = output_prefix+".core.phyml_tree_info.leaf_dists.group_dist.txt"

number_of_randomizations = 1000
total_jobs = 20

slurm_prefix = "#!/bin/bash\n#SBATCH -N 1#SBATCH -n 1\n#SBATCH -p sched_mit_chisholm,sched_mit_hill,newnodes\n#SBATCH --mem=1000\n#SBATCH --time=12:00:00\nmodule add engaging/R/3.1.1\ncd "+project_dir+"\n"

#############   FUNCTIONS   #############

def rand_groups(dict_in):
	key_list = []
	groups_dict = {}
	for key in dict_in:
		tup = (randint(0,10000),key)
		key_list.append(tup)
		try:
			groups_dict[dict_in[key]] += 1
		except:
			groups_dict[dict_in[key]] = 1
	key_list = sorted(key_list)
	#print(key_list)
	a = 0
	dict_out = {}
	group_dict_out = {}
	for group in groups_dict:
		out_string = ''
		count = groups_dict[group]
		for i in range(a,a+count):
			key = key_list[i][1]
			dict_out[key] = group
			if out_string == '':
				out_string = key
			else:
				out_string += ","+key
		group_dict_out[group] = out_string
		a = a+count
	# for group in group_dict_out:
	# 	print(group+"\t"+str(len(group_dict_out[group].split(","))))
	# print("\n")
	return group_dict_out


##################### MAIN #####################
##Make list of strains used in analysis
infile = open(project_dir+strain_list_file,"r")
master_strain_list = []
for line in infile:
	line = line.strip()
	line = line
	master_strain_list.append(line)
infile.close()

##Make dictionary of strain to populations/groups
pop_infile = open(project_dir+pop_infile_name,"r")
pop_dict = {}
pop_list = []
for line in pop_infile:
	line = line.strip().split("\t")
	strain = line[0]#+"_contigs"
	pop = line[1]
	pop_dict[strain] = pop
	pop_list.append(pop)
pop_infile.close()
pop_list = list(set(pop_list))

##Find maximum distance counted in each population
infile = open(input_dir+"max_distsum.txt","r")
max_dist_dict = {}
for line in infile:
	line = line.strip().split("\t")
	max_dist_dict[line[0]] = line[1]
infile.close()


iterations_per_job = (number_of_randomizations/total_jobs)+1
s = -1
for i in range(0,number_of_randomizations):
	rand_group_dict = rand_groups(pop_dict)
	group1_name = pop_list[0]
	group2_name = pop_list[1]
	group1 = rand_group_dict[group1_name]
	group2 = rand_group_dict[group2_name]
	group1_maxdist = max_dist_dict[group1_name]
	group2_maxdist = max_dist_dict[group2_name]

	if i == 0 or i%iterations_per_job == 0:
		s += 1
		mono_phy_sh = open(project_dir+"mono_phy_test."+str(s)+".sh","w")
		mono_phy_sh.write(slurm_prefix)
		mono_phy_sh.write("Rscript phybreak.monophy_test.R "+str(i)+" "+group1+" "+group2+"\n")
		mono_phy_sh.write("python phybreak.rand_group_leafdist.py "+str(i)+" "+group1+" "+group2+" "+group1_name+" "+group2_name+"\t"+group1_maxdist+"\t"+group2_maxdist+"\n")
		mono_phy_sh.write("rm "+input_dir+output_prefix+".rand_itt_"+str(i)+".group_monophy.txt\n")
	else:
		mono_phy_sh.write("Rscript phybreak.monophy_test.R "+str(i)+" "+group1+" "+group2+"\n")
		mono_phy_sh.write("python phybreak.rand_group_leafdist.py "+str(i)+" "+group1+" "+group2+" "+group1_name+" "+group2_name+"\t"+group1_maxdist+"\t"+group2_maxdist+"\n")
		mono_phy_sh.write("rm "+input_dir+output_prefix+".rand_itt_"+str(i)+".group_monophy.txt\n")
mono_phy_sh.close()

for i in range(0,s+1):
	command = "sbatch "+project_dir+"mono_phy_test."+str(i)+".sh"
	os.system(command)