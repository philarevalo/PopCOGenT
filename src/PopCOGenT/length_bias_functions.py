import numpy as np
from collections import Counter
from os import system, path, remove, makedirs
import random
import string
from Bio import SeqIO
from itertools import combinations, groupby


def align_and_calculate_length_bias(genome_1_file,
                                    genome_2_file,
                                    alignment_dir,
                                    mugsy_path,
                                    random_seed,
                                    keep_alignments):
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
    if not keep_alignments:
        remove(alignment_file)
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
    return mugsy_outname


def align_genomes(contig1,
                  contig2,
                  alignment_dir,
                  mugsy_path,
                  random_seed):
    
    random.seed(random_seed)
    # Assumes that files are named strain.contigextension.renamed.mugsy
    strain1 = '.'.join(path.basename(contig1).split('.')[0:-3])
    strain2 = '.'.join(path.basename(contig2).split('.')[0:-3])
    correct_name = '{strain1}_@_{strain2}.maf'.format(strain1 = strain1, strain2 = strain2) 
    final_name = alignment_dir+'/'+correct_name

    if not path.exists(final_name): # Only run the alignment if the file doesn't exist
        # make a temporary contig file due to parallelization issues with reading from the same filename
        out_id_1 = ''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for i in range(16))
        out_id_2 = ''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for i in range(16))


        system('cp {contig1} {alignment_dir}/{randomcontigname1}.tempcontig'.format(contig1=contig1, randomcontigname1=out_id_1, alignment_dir=alignment_dir))
        system('cp {contig2} {alignment_dir}/{randomcontigname2}.tempcontig'.format(contig2=contig2, randomcontigname2=out_id_2, alignment_dir=alignment_dir))

        # Aligning the genomes
        prefix = out_id_1 + out_id_2
        print('{mugsypath} --directory {align_directory} --prefix {prefix} {align_directory}/{randomcontigname1}.tempcontig {align_directory}/{randomcontigname2}.tempcontig'.format(mugsypath=mugsy_path,
                                                                                                                            align_directory=alignment_dir,
                                                                                                                            prefix = prefix,
                                                                                                                            randomcontigname1=out_id_1,
                                                                                                                            randomcontigname2 = out_id_2))
        system('{mugsypath} --directory {align_directory} --prefix {prefix} {align_directory}/{randomcontigname1}.tempcontig {align_directory}/{randomcontigname2}.tempcontig'.format(mugsypath=mugsy_path,
                                                                                                                            align_directory=alignment_dir,
                                                                                                                            prefix = prefix,
                                                                                                                            randomcontigname1=out_id_1,
                                                                                                                            randomcontigname2 = out_id_2))


        # Remove unneeded files
        remove('{align_directory}/{random_contig1}.tempcontig'.format(random_contig1=out_id_1, align_directory=alignment_dir))
        remove('{align_directory}/{random_contig2}.tempcontig'.format(random_contig2=out_id_2, align_directory=alignment_dir))
        remove('{align_directory}/{prefix}'.format(prefix=prefix, align_directory=alignment_dir))
        remove('{prefix}.mugsy.log'.format(prefix=prefix))

        system('mv {random_alignment_name} {correct_name}'.format(random_alignment_name=alignment_dir+'/'+prefix +'.maf',
                                                                  correct_name=alignment_dir+'/'+correct_name))
        
    return final_name

def calculate_length_bias(input_alignment,
                          genome_1_file,
                          genome_2_file,
                          output_file):


    g1size = sum([len(s) for s in SeqIO.parse(genome_1_file, 'fasta')])
    g2size = sum([len(s) for s in SeqIO.parse(genome_2_file, 'fasta')])

    if not path.exists(output_file): # only calculate the length bias if the file doesn't exist
        edge = get_transfer_measurement(input_alignment,
                                        g1size,
                                        g2size)

        with open(output_file, 'w') as outfile:
            outfile.write(edge + '\n')

def get_transfer_measurement(alignment,
                             g1size,
                             g2size,
                             min_block_size=0,
                             filtering_window=1000):

    # Initializes local variables
    filtered_blocks = []
    strain1, strain2 = alignment.split('/')[-1].split('_@_')
    strain2 = strain2.replace('.maf', '')
    all_blocks, prefilter_total_len = get_concatenated_alignment(alignment)

    # Filter alignment to split into subblocks at any point where there are at least 2 gaps
    for prefilter_s1, prefilter_s2 in all_blocks:
        filtered_blocks += filter_block(prefilter_s1, prefilter_s2)
    filtered_blocks = [block for block in filtered_blocks if len(block[0]) > min_block_size]
    s1temp, s2temp = zip(*filtered_blocks)

    # Assumes that each alignment block adds another divergence
    Concat_S1 = '1'.join(s1temp)
    Concat_S2 = '0'.join(s2temp)
    alignment_size = len(Concat_S1)
    init_div_count = naive_div_count(Concat_S1, Concat_S2)
    init_div = init_div_count * 1.0 / alignment_size

    # Second filtering step by divergence 
    final_filtered = []
    for s1, s2 in filtered_blocks:
        final_filtered += filter_block_by_divergence(s1, s2, init_div, winlen=filtering_window)
    filtered_blocks = [block for block in final_filtered if len(block[0]) > min_block_size]
    s1temp, s2temp = zip(*filtered_blocks)

    # Assumes that each alignment block adds another divergence
    Concat_S1 = '1'.join(s1temp)
    Concat_S2 = '0'.join(s2temp)
    alignment_size = len(Concat_S1)
    init_div_count = naive_div_count(Concat_S1, Concat_S2)
    init_div = (init_div_count * 1.0) / alignment_size

    if init_div > 0:
        initial = id_var_window_counts(Concat_S1, Concat_S2)
        initial_cumulative = get_cumulative_window_spectrum(initial, alignment_size)
        null_expect = single_param_null_model(np.arange(0, len(initial_cumulative)), init_div)
        observed_sum_sq_diff = np.sum(np.square(np.subtract(initial_cumulative, null_expect)))

        # Given a distribution of identical windows, bootstraps to find
        # length bias (SSD) confidence interval
        ssds = []
        for t in range(0, 200):
            initial_boot = np.random.choice(initial, len(initial), replace=True)
            initial_cumulative_boot = get_cumulative_window_spectrum(initial_boot, alignment_size)
            ssd_boot = np.sum(np.square(np.subtract(initial_cumulative_boot, null_expect)))
            ssds.append(ssd_boot)
        low_percentile = np.percentile(ssds, 0.5)
        high_percentile = np.percentile(ssds, 99.5)
    else:
        observed_sum_sq_diff = np.nan
        low_percentile = np.nan
        high_percentile = np.nan

    edge = '\t'.join([strain1,
                     strain2,
                     str(init_div),
                     str(alignment_size),
                     str(g1size),
                     str(g2size),
                     str(observed_sum_sq_diff),
                     str(low_percentile),
                     str(high_percentile)])
    return edge


def parse_alignment_file(alignment, min_block_size=1000, filtering_window=1000):

    # Initializes local variables
    filtered_blocks = []
    all_blocks, prefilter_total_len = get_concatenated_alignment(alignment)

    # Filter alignment to split into subblocks at any point where there are at least 2 gaps
    for prefilter_s1, prefilter_s2 in all_blocks:
        filtered_blocks += filter_block(prefilter_s1, prefilter_s2)
    filtered_blocks = [block for block in filtered_blocks if len(block[0]) > min_block_size]
    s1temp, s2temp = zip(*filtered_blocks)

    # Assumes that each alignment block adds another divergence
    Concat_S1 = '1'.join(s1temp)
    Concat_S2 = '0'.join(s2temp)
    alignment_size = len(Concat_S1)
    init_div_count = naive_div_count(Concat_S1, Concat_S2)
    init_div = init_div_count * 1.0 / alignment_size

    # Second filtering step by divergence 
    final_filtered = []
    for s1, s2 in filtered_blocks:
        final_filtered += filter_block_by_divergence(s1, s2, init_div, winlen=filtering_window)
    filtered_blocks = [block for block in final_filtered if len(block[0]) > min_block_size]
    s1temp, s2temp = zip(*filtered_blocks)

    # Assumes that each alignment block adds another divergence
    Concat_S1 = '1'.join(s1temp)
    Concat_S2 = '0'.join(s2temp)
    alignment_size = len(Concat_S1)
    init_div_count = naive_div_count(Concat_S1, Concat_S2)
    init_div = (init_div_count * 1.0) / alignment_size

    initial = id_var_window_counts(Concat_S1, Concat_S2)
    initial_cumulative = get_cumulative_window_spectrum(initial, alignment_size)

    return (initial_cumulative, init_div)


def get_cumulative_window_spectrum(idw, gs):
    '''
    Gets the X and Y coordinates of the identical window spectrum
    i.e., the fraction of the genome belonging to identical windows
    above a certain size
    '''

    obs_frac_counts = np.zeros(gs)
    norm = np.sum(idw)
    windows = Counter(idw)
    for wsize, count in windows.items():
        obs_frac_counts[wsize] = count * wsize * 1.0 / norm
    return 1.0 - np.cumsum(obs_frac_counts)


def get_concatenated_alignment(alignment):
    '''
    This creates a list of tuples that constitute a concatenated alignment.
    Every entry in the list is a tuple that corresponds to an alignment block.
    '''

    with open(alignment, 'r') as infile:
        '''
        Parser assumes a maf format where every alignment block begins with a
        statement of how many sequences are in that block, indicated by
        "mult=." Also assumes that the order of sequences in each block is
        the same.
        '''
        seqs = []
        total_len = 0
        for lines in infile:
            if 'mult=2' in lines:
                seq_line_1 = next(infile)
                block_1 = seq_line_1.split()[-1].strip()
                total_len += len(block_1)
                seq_line_2 = next(infile)
                block_2 = seq_line_2.split()[-1].replace('\n', '')
                seqs.append((block_1, block_2))
    return seqs, total_len


def id_var_window_counts(sequence_1, sequence_2):
    '''
    This method takes two aligned sequences (strings) and returns the
    lengths of all identical windows between those sequences.
    '''
    if sequence_1 == sequence_2:
        id_seqs = [len(sequence_1)]
    else:
        a1 = np.array(list(sequence_1))
        a2 = np.array(list(sequence_2))
        mutated_positions = np.where(a1 != a2)[0]
        id_seqs = -1 + np.ediff1d(mutated_positions,
                                  to_begin=mutated_positions[0] + 1,
                                  to_end=len(sequence_1) - mutated_positions[-1])
    return id_seqs


def naive_div_count(sequence_1, sequence_2):
    '''
    Given two aligned strings, returns the number of differences between them.
    '''
    if sequence_1 == sequence_2:
        return 0
    a1 = np.array(list(sequence_1))
    a2 = np.array(list(sequence_2))
    return len(np.where(a1 != a2)[0])


def filter_block(sequence_1, sequence_2):
    removal_positions = filter_string(sequence_1)
    removal_positions += filter_string(sequence_2)
    return [block for block in get_filtered_subblocks(sequence_1, sequence_2, removal_positions) if block[0] != '']


def filter_string(S):
    groups = groupby(S)
    result = [(label, sum(1 for _ in group)) for label, group in groups]
    begin = 0
    filter_intervals = []
    for base, count in result:
        end = begin + count
        if base == '-' and count >= 2:
            filter_intervals.append((end, begin))
        begin += count
    return(filter_intervals)


def filter_block_by_divergence(sequence_1, sequence_2, init_div, winlen=1000):
    '''
    Filters two sequences from an alignment block to remove regions
    that are significantly more diverged than expected
    '''

    if sequence_1 == sequence_2:
        return [(sequence_1, sequence_2)]
    else:
        removal_positions = []
        begin = 0
        for end in range(winlen, len(sequence_1), winlen):
            d = naive_div_count(sequence_1[begin:end], sequence_2[begin:end])
            if d / winlen >= 10 * init_div:
                removal_positions.append((end, begin))
            begin = end

        if begin < len(sequence_1):
            d = naive_div_count(sequence_1[begin:], sequence_2[begin:])
            if d / (len(sequence_1) - begin) >= 10 * init_div:
                removal_positions.append((len(sequence_1), begin))

        return [block for block in get_filtered_subblocks(sequence_1, sequence_2, removal_positions) if block[0] != '']


def get_filtered_subblocks(sequence_1, sequence_2, positions_to_remove):
    '''
    Helper method that splits a string into substrings when given a list of
    start and end positions to remove
    '''
    if positions_to_remove == []:
        return [(sequence_1, sequence_2)]
    else:
        final_blocks = []
        initial_start = 0
        if len(positions_to_remove) > 1:
            merged = merge_intervals(sorted(positions_to_remove, reverse=True))
        else:
            merged = positions_to_remove
        ends, starts = zip(*sorted(merged, key=lambda x: x[1]))
        for i, start_of_deleted_region in enumerate(starts):
            end_of_deleted_region = ends[i]
            subsequence_1 = sequence_1[initial_start:start_of_deleted_region]
            subsequence_2 = sequence_2[initial_start:start_of_deleted_region]
            initial_start = end_of_deleted_region
            final_blocks.append((subsequence_1, subsequence_2))
        final_blocks.append((sequence_1[initial_start:], sequence_2[initial_start:]))
        return final_blocks


def single_param_null_model(w, div):
    '''
    The simple single parameter null model that describes
    the window spectrum under an assumption of only mutation
    and no transfer
    '''
    return np.exp(-div * w) * (div * w + 1)


def merge_intervals(intervals):
    all_intervals = []
    for j, interval in enumerate(intervals):
        end, start = interval
        if j == 0:
            current_interval = interval
        else:
            if intervals[j-1][1] <= end:
                current_interval = (current_interval[0], start)
            else:
                all_intervals.append(current_interval)
                current_interval = interval
    if len(all_intervals) > 0:
        if all_intervals[-1] != current_interval:
            all_intervals.append(current_interval)
    else:
        all_intervals.append(current_interval)
    return all_intervals

