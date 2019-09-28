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

#############   MAIN   #############

##Make dictionary of strain to populations/groups
pop_infile = open(project_dir+pop_infile_name,"r")
pop_dict = {}
pop_list = []
for line in pop_infile:
    line = line.strip().split("\t")
    strain = line[0]
    pop = line[1]
    if strain != 'Strain':
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
infile = open(input_dir+leaf_dist_file,"r")
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
                if pop1 == focus_population:
                    pop1 = "focus"
                else:
                    pop1 = "other"
                pop2 = pop_dict[iso2]
                if pop2 == focus_population:
                    pop2 = "focus"
                else:
                    pop2 = "other"
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

##Collect tree locations in the MSA
infile = open(input_dir+tree_info_file,"r")
tree_loc_info = {}
for line in infile:
    line = line.strip().split("\t")
    tree_no = line[0]
    msa_start = int(line[1])
    msa_stop = int(line[2])
    tree_loc_info[tree_no] = (msa_start,msa_stop)
infile.close()


##Output the fraction of total branch length within a popultion for each tree
percentile_dict = {}
physplit_trees = []
focus_other_list = ['other','focus']
with open(phy_split_dir+summary_file_name,"w") as outfile:
    outfile.write("mid_point_of_window\ttree_no\tmonophy\tfocus\n")
    for j in range(1,last_tree_no+1): # in dist_dict:
        tree_no = str(j)

        mid_point_of_window = (tree_loc_info[tree_no][0]+tree_loc_info[tree_no][1])/2
        
        dist_list = dist_dict[tree_no]['focus\tfocus']
        dist_sum = str(numpy.sum(dist_list)/all_pairwise_sum[tree_no]) # intrapopulation distance fraction for tree of interest
        outfile.write(str(mid_point_of_window) + "\t" + tree_no + "\t" + mono_phy_dict[tree_no] + "\t" + dist_sum + "\n")