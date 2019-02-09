import numpy as np
from collections import Counter
from joblib import Parallel, delayed
import itertools
import argparse
from Bio import SeqIO
from length_bias_functions import *


def main():
    parser = argparse.ArgumentParser(
        description=('Calculate length bias as the sum of squared differences between aligned genomes.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i',
                        '--input_alignment',
                        help='Input alignment file.')
    parser.add_argument('-g1',
                        '--genome_1_file',
                        help='Path to genome 1.')
    parser.add_argument('-g2',
                        '--genome_2_file',
                        help='Path to genome 2')
    parser.add_argument('-o',
                        '--output_file',
                        help='Output file name')
    parser.set_defaults(simulate_transfer=False)

    args = parser.parse_args()

    g1size = sum([len(s) for s in SeqIO.parse(args.genome_1_file, 'fasta')])
    g2size = sum([len(s) for s in SeqIO.parse(args.genome_2_file, 'fasta')])

    edge = get_transfer_measurement(args.input_alignment,
                                    g1size,
                                    g2size,
                                    args.simulate_transfer)

    with open(args.output_file, 'w') as outfile:
        outfile.write(edge + '\n')


def get_transfer_measurement(alignment,
                             g1size,
                             g2size,
                             min_block_size=1000,
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


if __name__ == '__main__':
    main()
