from Bio import SeqIO
import glob
import pandas as pd
from itertools import combinations
import os
from length_bias_functions import *


configfile:
    "PopCOGenT_config.yml"

genome_directory = config["genome_directory"]
genome_extension = config["genome_extension"]
alignment_dir = config["alignment_dir"]
base_name = config["base_name"]
output_dir = config["output_directory"]
mugsy_path = config["mugsy_path"]
clonal_cutoff = config["clonal_cutoff"]

strains = [os.path.splitext(os.path.basename(g))[0] for g in glob.glob('%s/*%s'%(genome_directory, genome_extension))]
combos = combinations(strains, 2)
strains1, strains2 = zip(*combos)
seeds = {(g1, g2): random.randint(1, int(1e9)) for g1, g2 in combinations([g for g in glob.glob('%s/*%s'%(genome_directory, genome_extension))], 2)}

rule target:    
    input:
        '%s/%s_%s.txt.cluster.tab.txt'%(output_dir, base_name, clonal_cutoff)

rule cluster:
    input:
        "%s/%s.length_bias.txt"%(output_dir, base_name)
    output:
        '%s/%s_%s.txt.cluster.tab.txt'%(output_dir, base_name, clonal_cutoff)

    shell:
        "python cluster.py --base_name {config[base_name]} --length_bias_file {input} --clonal_cutoff {config[clonal_cutoff]} --infomap_path {config[infomap_path]} --output_directory {config[output_directory]}"

rule concatenate_length_bias_files:
    input:
        expand(alignment_dir + "{g1}_@_{g2}.length_bias.txt", zip, g1=strains1, g2=strains2)
    output:
        "%s/%s.length_bias.txt"%(output_dir, base_name)
    run:

        header = ['Strain 1',
         'Strain 2',
         'Initial divergence',
         'Alignment size',
         'Genome 1 size',
         'Genome 2 size',
         'Observed SSD',
         'SSD 95 CI low',
         'SSD 95 CI high']
        rows = [open(f).read().strip().split() for f in input]
        df = pd.DataFrame(rows, columns=header)
        df.to_csv(str(output), sep='\t', index=False)

rule make_length_bias_file:
    input:
        g1 = genome_directory + "/{g1}" + genome_extension,
        g2 = genome_directory + "/{g2}" + genome_extension
    output: 
        alignment_dir + "{g1}_@_{g2}.length_bias.txt"
    run:
        seed = seeds[(input['g1'], input['g2'])]
        g1 = rename_for_mugsy(input['g1'])
        g2 = rename_for_mugsy(input['g2'])
        aln = align_genomes(g1, g2, alignment_dir, mugsy_path, seed)

        calculate_length_bias(aln, g1, g2, str(output))