import os
import sys
###########   USER DEFINED VARIABLES   #############
project_dir = ""
contig_extension = ""
strain_list_file = ""
output_prefix = ""

with open("phybreak_parameters.txt","r") as parameter_file:
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
input_dir = project_dir+"align/"
alignment_dir = input_dir+"alignment_blocks/"
phy_split_dir = input_dir+"phy_split/"

MSA_name = output_prefix+".core.fasta"
phy_prefix = output_prefix+".core"
block_loc = output_prefix+".block_location.txt"
snp_loc_file = output_prefix+".core.SNPloc.txt"

window_size = 100 #number of SNPs to include per tree
overlap = 1 #number of SNPs to overlap between trees

total_jobs = 50 #number of jobs total to generate

slurm_prefix = "#!/bin/bash\n#SBATCH -N 1#SBATCH -n 1\n#SBATCH -p sched_mit_chisholm,sched_mit_hill,newnodes\n#SBATCH --mem=2000\n#SBATCH --time=48:00:00\ncd "+project_dir+"\n"

#phyML sh file variables
phyML_loc = "PhyML "
phyML_properties = "-q -m JC69 -f e -c 2 -a 0.022" #-m JC69 -f e -c $Ncat -a $alpha 
#############   FUNCTIONS   #############

def check_range_overlap(range_start,range_stop,range_list):
	outval = True
	for break_point in range_list:
		if break_point >= range_start and break_point <= range_stop:
			outval = False
			break
	return outval

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

#############   MAIN   #############

if os.path.isdir(phy_split_dir) == False:
	os.makedirs(phy_split_dir)

#store MSA in dictionary
with open(input_dir+MSA_name,"r") as msa:
	seq_dict = {}
	head = ""
	for line in msa:
	    line = line.strip()
	    if line[0] == ">":
	        head = line[1:len(line)]
	    else:
	        try:
	            seq_dict[head] += line
	        except:
	            seq_dict[head] = line
msa_len = len(seq_dict[head])

print("Done storing MSA")

#find locations of blocks in alignment
with open(input_dir + block_loc,"r") as block_file:
	break_list = []
	a = 0
	for line in block_file:
		a += 1
		if a > 1:
			line = line.strip().split("\t")
			break_list.append(int(line[3]))
print("Done storing LCA locations")

#store locations of SNPs in dictionary
with open(input_dir+snp_loc_file,"r") as SNPfile:
	snp_dict = {}
	for line in SNPfile:
		line = line.strip().split("\t")
		snp_dict[int(line[0])] = int(line[1])

total_snps = len(snp_dict)
print("Done storing SNP locations")

####count number of trees we will make
print("Searching for number of total trees to make")
tree_total = 0
end = 0
while end < total_snps-window_size:
	strt = snp_dict[end]
	stp = snp_dict[end+window_size]
	window_length = stp-strt
	if check_range_overlap(strt,stp,break_list) == True:
		tree_total += 1
		end += overlap
	else:
		end += 1
		#print(end)
trees_per_file = int(float(tree_total)/float(total_jobs))
print("Total trees: "+str(tree_total))

####move across MSA and write phylip phy files
end = 0
tree_count = 0
phy_num = 0
phy_count = {} #dictionary of number of trees per PHYLIP formatted file
while end < total_snps-window_size:
	strt = snp_dict[end]
	stp = snp_dict[end+window_size]
	window_length = stp-strt
	if check_range_overlap(strt,stp,break_list) == True:
		subset = {}
		for head in seq_dict:
			seq = seq_dict[head][strt:stp]
			subset[head] = seq
		phylip_seq = fasta_2_phylip(subset,window_length)
		end += overlap
		tree_count += 1

		# writes out alignment windows individually
		with open(phy_split_dir + phy_prefix + "." + str(tree_count) + "window.phy", "w") as window_out:
			window_out.write(phylip_seq)

		if tree_count == 1:
			outfile = open(phy_split_dir+phy_prefix+"."+str(phy_num)+".phy","w")
			outfile.write(phylip_seq)

			treeloc = open(phy_split_dir+phy_prefix+"."+str(phy_num)+".treeloc.txt","w")
			treeloc.write(str(tree_count) +"\t"+ str(strt) +"\t"+ str(stp)+"\n")
			
			phy_count[phy_num] = 1
		elif tree_count % trees_per_file == 0:
			outfile.write(phylip_seq)
			phy_count[phy_num] += 1
			
			phy_num += 1
			phy_count[phy_num] = 0
			outfile = open(phy_split_dir+phy_prefix+"."+str(phy_num)+".phy","w")
			
			treeloc.write(str(tree_count) +"\t"+ str(strt) +"\t"+ str(stp)+"\n")
			treeloc = open(phy_split_dir+phy_prefix+"."+str(phy_num)+".treeloc.txt","w")

		else:
			outfile.write(phylip_seq)
			treeloc.write(str(tree_count) +"\t"+ str(strt) +"\t"+ str(stp)+"\n")
			phy_count[phy_num] += 1

	else:
		end += 1
outfile.close()
treeloc.close()

#write SLURM files for phyML job sbatch
for num in phy_count:
	count = str(phy_count[num])
	with open(project_dir+phy_prefix+"."+str(num)+".phyML.sh","w") as phyml_sh:
		phy_line = phyML_loc+"-i "+phy_split_dir+phy_prefix+"."+str(num)+".phy -n " +count+ " "
		phy_line += phyML_properties +" > "+phy_split_dir+phy_prefix+"."+str(num)+".phy_phyml_stat.txt\n"
		phyml_sh.write(slurm_prefix+phy_line)
	command = "sbatch "+project_dir+phy_prefix+"."+str(num)+".phyML.sh"
	os.system(command)

print(str(tree_count)+" trees written to "+str(len(phy_count))+" files")
