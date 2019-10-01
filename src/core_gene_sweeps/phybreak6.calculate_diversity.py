import pandas as pd
from Bio import SeqIO
import glob
from itertools import combinations, product
from collections import defaultdict
import numpy as np
import os

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

input_dir = project_dir+"align/"
phy_split_dir = input_dir+"phy_split/"


def count_divs(s1, s2):
    '''
    Given two strings counts the number of differences between them
    '''
    d = 0
    for b1, b2 in zip(s1, s2):
        if b1 != b2:
            d += 1
    return d * 1.0 / len(s1)


def calc_pop_div(pop, s_dict):
    '''
    Given a list of strains and a dictionary of sequences, calculates pi,
    the average nucleotide diversity between pairs of strains
    '''
    
    # pop is just a list of strains in the population of interest
    # s dict is a dictionary that maps strains/genomes to sequences
    aves = []
    for strain1, strain2 in combinations(pop, 2):
        if strain1 in s_dict.keys() and strain2 in s_dict.keys():
            aves.append(count_divs(s_dict[strain1], s_dict[strain2]))
            L = len(s_dict[strain1])
    return(np.average(aves), L)

def clean_id(clone):
    c = clone.strip().split(',')[0]
    return c

os.system('mkdir -p output')
df = pd.read_table(pop_infile_name, dtype=str)

# Loops over all alignment blocks and puts sequences into a dictionary
seqs = {}
for s in SeqIO.parse('%s/%s.core.fasta'%(input_dir, output_prefix), 'fasta'):
    seqs[s.id] = str(s.seq)

print('Calculating genome-wide diversity.')
# Calculate pi for each population
new_rows = []
for pop, popdf in df.groupby('Cluster_ID'):
    strains = [clean_id(clone) for clone in set(popdf.Clonal_complex)]
    pi, length = calc_pop_div(strains, seqs) 
    new_rows.append([pop, pi])
population_pi = pd.DataFrame(new_rows, columns=['Cluster_ID', 'Pi'])
population_pi.to_csv('output/%s_pop_pi.csv'%output_prefix, index=False)
print('Done calculating genome-wide diversity.')

print('Calculating diversity of individual alignments')
# Split alignments of all 100 SNP trees into individual trees
tree_no = 1
new_lines = []
new_rows = []
for line in open('%s/%s.phy'%(phy_split_dir, output_prefix)):
    if line == '\n':
        if len(new_lines) > 0:
            all_seqs = {l.split()[0].strip(): l.split()[-1].strip() for l in new_lines}
            for pop, popdf in df.groupby('Cluster_ID'):
                strains = [clean_id(clone) for clone in set(popdf.Clonal_complex)]
                # Diversity within a population
                intra_pop_div, length = calc_pop_div(strains, all_seqs)
                new_rows.append([tree_no, pop, intra_pop_div, length])
            if tree_no%100 == 0:
                print('Finished calculation for %s alignments'%tree_no)
            tree_no += 1
        new_lines = []
    else:
        new_lines.append(line)

print('Done calculating individual alignment diversity')
intra_df = pd.DataFrame(new_rows,columns=['tree_no', 'Cluster_ID', 'block_pi', 'Length'])
intra_df.to_csv('output/%s_block_pi.csv'%output_prefix, index=False)
