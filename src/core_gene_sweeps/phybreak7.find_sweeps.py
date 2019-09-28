from Bio import SeqIO
import glob
import pandas as pd
from itertools import combinations, product
import numpy as np
from scipy import stats 
import os


def concatenate_windows(sweepdf, length_cutoff=1000):
    i = 0
    final_df = pd.DataFrame()
    current_end = None
    blocks = []
    for index, row in sweepdf.iterrows():
        if not current_end:
            current_end = row['End']
            blocks.append(index)
        elif row['Start'] <= current_end:
            blocks.append(index)

        else: 
            if len(blocks) > 0:
                if sweepdf.loc[max(blocks), 'End'] - sweepdf.loc[min(blocks), 'Start'] >= length_cutoff:
                    block_begin = min(blocks)
                    block_end = max(blocks)
                    final_df.loc[i, 'Start'] = sweepdf.loc[min(blocks), 'Start']
                    final_df.loc[i, 'End'] = sweepdf.loc[max(blocks), 'End']
                    final_df.loc[i, 'Start tree'] = min(blocks)
                    final_df.loc[i, 'End tree'] = max(blocks)
                    i += 1
            blocks = []
        current_end = row['End']
        
    if len(blocks) > 0:
        if sweepdf.loc[max(blocks), 'End'] - sweepdf.loc[min(blocks), 'Start'] >= length_cutoff:
            block_begin = min(blocks)
            block_end = max(blocks)
            final_df.loc[i, 'Start'] = sweepdf.loc[min(blocks), 'Start']
            final_df.loc[i, 'End'] = sweepdf.loc[max(blocks), 'End']
            final_df.loc[i, 'Start tree'] = min(blocks)
            final_df.loc[i, 'End tree'] = max(blocks)
    if len(final_df) > 0:
        final_df['Midpoint'] = final_df.Start + ((final_df.End - final_df.Start) // 2)
    else:
        print('No sweeps found')
    return final_df

def calc_all_divs(concatenated_positions,
                  population,
                  pi,
                  alpha):
    finaldf = pd.DataFrame(index=concatenated_positions.index, columns=['Start',
                                                                        'End',
                                                                        'Start tree',
                                                                        'End tree',
                                                                        'pop_pi'])
    for i in concatenated_positions.index:
        start = int(concatenated_positions.loc[i, 'Start'])
        end = int(concatenated_positions.loc[i, 'End'])
        seqs = {}
        for strain, s in total_alignment.items():
            seqs[strain] = s[start: end]
        
        pop_div = calc_pop_div(population, seqs)

        finaldf.loc[i, ] = [start,
                            end,
                            concatenated_positions.loc[i, 'Start tree'],
                            concatenated_positions.loc[i, 'End tree'],
                            pop_div]
    finaldf['Length'] = finaldf.End - finaldf.Start
    finaldf['Midpoint'] = finaldf.Start + ((finaldf.End - finaldf.Start) // 2)
    calculate_ci(finaldf, alpha, pi)
    return finaldf
    
def count_divs(s1, s2):
    d = 0
    for b1, b2 in zip(s1, s2):
        if b1 != b2:
            d += 1
    return d * 1.0 / len(s1)

def calc_pop_div(pop, s_dict):
    aves = []
    for strain1, strain2 in combinations(pop, 2):
        if strain1 in s_dict.keys() and strain2 in s_dict.keys():
            aves.append(count_divs(s_dict[strain1], s_dict[strain2]))
    return(np.average(aves))

def calculate_ci(df, alpha, p_pop):
    low, high = stats.binom.interval(alpha, list(df.Length), p_pop)
    df['ci_low'] = low / df.Length
    df['ci_high'] = high / df.Length

def passes_phylo_criteria(full_df, cutoff):
    return (full_df.focus < cutoff) & (full_df.monophy == 1) & ((full_df.block_pi < full_df.ci_low) | (full_df.block_pi==0))

def find_sweeps(phybreak_df,
                pi_info_df,
                pop_pi, alpha,
                pop_name,
                pop_list,
                threshold):
    full_df = pd.merge(pi_info_df, phybreak_df)
    full_df.index = full_df.block
    calculate_ci(full_df, alpha, pop_pi)
    # Finds cutoff based on fraction of branch length within the focus population
    cutoff = np.percentile(full_df.focus, threshold)
    trees_passing_phy_criteria = full_df[passes_phylo_criteria(full_df, cutoff)]
    concat = concatenate_windows(trees_passing_phy_criteria, length_cutoff=500)
    concat = calc_all_divs(concat, pop_list, pop_pi, alpha)
    for col in concat:
        concat[col] = concat[col].astype(float)
    return concat



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


os.system('mkdir -p output')


intra_pi = pd.read_csv('output/%s_block_pi.csv'%output_prefix, dtype=str)

pop_pi = pd.read_csv('output/%s_pop_pi.csv'%output_prefix, index_col=0, dtype=str)
pop_pi.index = pop_pi.index.map(str)

phybreak_result = pd.read_csv('%s/phybreak_result_%s.txt'%(phy_split_dir, focus_population), delimiter='\t')

treeloc = pd.read_csv('%s/%s.treeloc.txt'%(phy_split_dir, output_prefix), delimiter='\t', header=None, dtype=str)
treeloc.columns = ['tree_no', 'Start', 'End']

msa_file = '%s/%s.core.fasta'%(input_dir, output_prefix)
total_alignment = {s.id : str(s.seq) for s in SeqIO.parse(msa_file, 'fasta')}

pops = pd.read_table(pop_infile_name, dtype=str)
pop_list = list(pops[pops.Cluster_ID == focus_population]['Strain'])

seqs = {}
for s in SeqIO.parse(msa_file, 'fasta'):
    seqs[s.id] = str(s.seq)

p0 = float(pop_pi.loc[focus_population, 'Pi'])
div_df = intra_pi[intra_pi.Cluster_ID == focus_population].merge(treeloc)
for col in div_df.columns:
    div_df[col] = div_df[col].astype(float)
div_df['block'] = div_df['tree_no']
div_df.index = div_df['tree_no']
div_df['Midpoint'] = div_df.Start + (div_df.End - div_df.Start) // 2

div_df.index.name = ''
concat_pop = find_sweeps(phybreak_result, div_df, p0, 0.95, focus_population, pop_list, percentile_threshold)
final = concat_pop[(concat_pop.pop_pi < concat_pop.ci_low) | (concat_pop.pop_pi == 0)]
final.to_csv('output/%s.%s.core_sweeps.csv'%(output_prefix, focus_population))