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
msa = open(input_dir+MSA_name,"r")
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
msa.close()
msa_len = len(seq_dict[head])
print("Done storing MSA")

#find locations of blocks in alignment
block_file = open(input_dir+block_loc_filename,"r")
break_list = []
a = 0
for line in block_file:
	a += 1
	if a > 1:
		line = line.strip().split("\t")
		break_list.append(int(line[3]))
print("Done storing LCB locations")

#store locations of SNPs in dictionary
SNPfile = open(input_dir+snp_loc_filename,"r")
snp_dict = {}
for line in SNPfile:
	line = line.strip().split("\t")
	snp_dict[int(line[0])] = int(line[1])
SNPfile.close()
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
print("Total trees: "+str(tree_total))

#### move across MSA and write phylip format .phy files
end = 0
tree_count = 0
phy_count = 0 #dictionary of number of trees per PHYLIP formatted file
phylip_outfile = open(phy_split_dir+phy_prefix+".phy","w")
treeloc = open(phy_split_dir+treeloc_filename,"w")
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
		phylip_outfile.write(phylip_seq)

		treeloc.write(str(tree_count) +"\t"+ str(strt) +"\t"+ str(stp)+"\n")
		
		phy_count += 1
	else:
		end += 1
phylip_outfile.close()
treeloc.close()

#run phyML to generate trees
phy_command = phyML_loc+" -i "+phy_split_dir+phy_prefix+".phy -n " +str(phy_count)+ " "
phy_command += phyML_properties +" > "+phy_split_dir+phy_prefix+".phy_phyml_stat.txt\n"
print(phy_command)
os.system(phy_command)

