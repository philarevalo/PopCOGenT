import argparse
from os import system, path, remove
import random
import string
from Bio import SeqIO
from length_bias_functions import *
from joblib import Parallel, delayed
import glob
from itertools import combinations


def main():
    run_on_single_machine(10,
                         '../../test/',
                         '.fasta',
                         '../../test/',
                         '~/apps/mugsy_trunk/mugsy')



def run_on_single_machine(threads,
                          genome_directory,
                          contig_extension,
                          alignment_dir,
                          mugsy_path):
    pairs_and_seeds = [(g1, g2, random.randint(1, int(1e9))) for g1, g2 in combinations(glob.glob(genome_directory + '*' + contig_extension), 2)]
    Parallel(n_jobs=threads)(delayed(align_and_calculate_length_bias)(g1, g2, alignment_dir, mugsy_path, seed) for g1, g2, seed in pairs_and_seeds)

def align_and_calculate_length_bias(genome_1_file,
                                    genome_2_file,
                                    alignment_dir,
                                    mugsy_path,
                                    random_seed):
    alignment_file = align_genomes(genome_1_file,
                                   genome_2_file,
                                   alignment_dir,
                                   mugsy_path,
                                   random_seed)

    length_bias_file = alignment_file + '.length_bias.txt'
    calculate_length_bias(alignment_file,
                          genome_1_file,
                          genome_2_file,
                          length_bias_file)
    return length_bias_file



def rename_for_mugsy(genome):

    # Assumes the strain name is everythign except the extension
    strain_name = '.'.join(path.basename(genome).split('.')[0:-1])

    # We want to remove all periods and colons from sequence input so that mugsy doesn't break
    mugsy_outname = genome + '.renamed.mugsy'
    
    # Removes all bad characters
    mugsy_name = strain_name.translate(({ord(c): '_' for c in """ !@#$%^&*()[]{};:,./<>?\|`"'~-=+"""}))
    
    mugsy_s = []
    for i, s in enumerate(SeqIO.parse(genome, 'fasta')):
        s.description = ''
        s.id = '{id}_{contig_num}'.format(id=mugsy_name, contig_num=str(i))
        mugsy_s.append(s)
    SeqIO.write(mugsy_s, mugsy_outname, 'fasta')


def align_genomes(contig1,
                  contig2,
                  alignment_dir,
                  mugsy_path,
                  random_seed):
    
    random.seed(random_seed)
    # Assumes that files are named strain.contigextension.renamed.mugsy
    strain1 = path.basename.split('.')[0:-3]
    strain2 = path.basename.split('.')[0:-3]
    correct_name = '{strain1}_@_{strain2}.maf'.format(strain1 = strain1, strain2 = strain2) 

    # make a temporary contig file due to parallelization issues with reading from the same filename
    out_id_1 = ''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for i in range(16))
    out_id_2 = ''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for i in range(16))


    system('cp {contig1} {randomcontigname1}'.format(contig1=contig1, randomcontigname1=out_id_1))
    system('cp {contig2} {randomcontigname2}'.format(contig2=contig2, randomcontigname2=out_id_2))

    # Aligning the genomes
    prefix = out_id_1 + out_id_2
    system('{mugsypath} --directory {align_directory} --prefix {prefix} {randomcontigname1} {randomcontigname2}'.format(mugsypath=mugsy_path,
                                                                                                                        align_directory=alignment_dir,
                                                                                                                        prefix = prefix,
                                                                                                                        randomcontigname1=out_id_1,
                                                                                                                        randomcontigname2 = out_id_2))

    # Remove unneeded files
    remove('{random_contig1}'.format(random_contig1=out_id_1))
    remove('{random_contig2}'.format(random_contig2=out_id_2))
    remove('{prefix}.mugsy.log'.format(prefix=prefix))

    system('mv {random_alignment_name} {correct_name}'.format(random_alignment_name=alignment_dir+'/'+prefix +'.maf',
                                                              correct_name=alignment_dir+'/'+correct_name))
    return alignment_dir+'/'+correct_name

def calculate_length_bias(input_alignment,
                          genome_1_file,
                          genome_2_file,
                          output_file):


    g1size = sum([len(s) for s in SeqIO.parse(genome_1_file, 'fasta')])
    g2size = sum([len(s) for s in SeqIO.parse(genome_2_file, 'fasta')])

    edge = get_transfer_measurement(input_alignment,
                                    g1size,
                                    g2size)

    with open(output_file, 'w') as outfile:
        outfile.write(edge + '\n')