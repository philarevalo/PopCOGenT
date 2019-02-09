##This script consolodates all of the information from MSA_subset_and_writePHY.py

###########   USER DEFINED VARIABLES   #############
project_dir = ""
contig_extension = ""
strain_list_file = ""
output_prefix = ""
window_size = 0 #number of SNPs to include per tree

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
        elif line[0] == "window_size":
            window_size = int(line[1].split(" #")[0])
        # elif line[0] == "ref_iso":
        #     ref_iso = line[1].split(" #")[0]
        # elif line[0] == "ref_contig":
        #     ref_contig = line[1].split(" #")[0]
        # elif line[0] == "contig_dir":
        #     contig_dir = line[1].split(" #")[0]
        # elif line[0] == "pop_infile_name":
        #     pop_infile_name = line[1].split(" #")[0]
        # elif line[0] == "len_block_threshold":
        #     len_block_threshold = int(line[1].split(" #")[0])
        # elif line[0] == "gap_prop_thresh":
        #     gap_prop_thresh = float(line[1].split(" #")[0])
parameter_file.close()

overlap = 1 #number of SNPs to overlap between trees
total_jobs = 50 #number of jobs generated in last step

input_dir = project_dir+"align/"
alignment_dir = input_dir+"alignment_blocks/"
phy_split_dir = input_dir+"phy_split/"
tree_dir = input_dir+"trees/"

MSA_name = output_prefix+".core.fasta"
phy_prefix = output_prefix+".core"


#############   FUNCTIONS   #############

#############   MAIN   #############
import os
import sys
if os.path.isdir(tree_dir) == False:
    os.makedirs(tree_dir)


ML_dict = {}
for i in range(0,total_jobs+1):
    subseq = phy_prefix+"."+str(i)
    #print(subseq)

    ##Build dictionary of PhyML  likelihood values
    lkfile = open(phy_split_dir+subseq+".phy_phyml_stat.txt","r")
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
                        ML_dict[subseq][tree_no] = lk +"\t"+ nonsnp #+"\t"
                    except:
                        ML_dict[subseq] = {}
                        ML_dict[subseq][tree_no] = lk +"\t"+ nonsnp #+"\t"
                    tree_no = ""
                    lk = ""
    lkfile.close()

    ##Add location of tree in MSA to dictionary
    treeloc_file = open(phy_split_dir+subseq+".treeloc.txt","r")
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
outfile = open(input_dir+phy_prefix+".phyml_tree_info.txt","w")
for i in range(0,total_jobs+1):
    subseq = phy_prefix+"."+str(i)
    for tree_no in range(1,len(ML_dict[subseq])+1):
        outfile.write(ML_dict[subseq][str(tree_no)]+"\n")
        #write tree to file
        tree_num = ML_dict[subseq][str(tree_no)].split("\t")[0]
        tree_string = ML_dict[subseq][str(tree_no)].split("\t")[5]
        # tree_outfile = open(tree_dir+str(tree_num)+".nwk","w")
        # tree_outfile.write(tree_string+"\n")
        # tree_outfile.close()
outfile.close()
