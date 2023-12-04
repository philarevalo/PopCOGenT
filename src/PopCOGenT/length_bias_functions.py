import os
import random
import string
from collections import Counter
from itertools import groupby
from typing import Iterable, List, Tuple

import numpy as np
from Bio import SeqIO


def align_and_calculate_length_bias(
    genome_1_file,
    genome_2_file,
    alignment_dir,
    mugsy_path,
    random_seed,
    keep_alignments,
):
    alignment_file = align_genomes(
        genome_1_file, genome_2_file, alignment_dir, random_seed
    )

    length_bias_file = alignment_file + ".length_bias.txt"
    calculate_length_bias(
        alignment_file, genome_1_file, genome_2_file, length_bias_file
    )
    if not keep_alignments:
        os.remove(alignment_file)
    return length_bias_file


def rename_for_mugsy(genome, out_genome=None, strain_name=None):
    # We want to remove all periods and colons from sequence input so that mugsy doesn't break
    mugsy_outname = out_genome or (genome + ".renamed.mugsy")
    # Assumes the strain name is everythign except the extension
    strain_name = strain_name or os.path.basename(genome).rsplit(".", 1)[0]

    # Removes all bad characters
    mugsy_name = strain_name.translate(
        {ord(c): f"_ord{ord(c)}_" for c in """ !@#$%^&*()[]{};:,./<>?\|`"'~-=+"""}
    )

    def _yield():
        for i, s in enumerate(SeqIO.parse(genome, "fasta")):
            s.description = ""
            s.id = f"{mugsy_name}_{i}"
            yield s

    SeqIO.write(_yield(), mugsy_outname, "fasta-2line")
    return mugsy_outname


def align_genomes(contig1: str, contig2: str, alignment_dir: str, random_seed: int):
    # Assumes that files are named strain.contigextension.renamed.mugsy
    strain1 = os.path.basename(contig1).rsplit(".", 3)[0]
    strain2 = os.path.basename(contig2).rsplit(".", 3)[0]
    correct_name = f"{strain1}_@_{strain2}.maf"
    final_name = f"{alignment_dir}/{correct_name}"

    if not os.path.exists(
        final_name
    ):  # Only run the alignment if the file doesn't exist
        random.seed(random_seed)

        # make a temporary contig file due to parallelization issues with reading from the same filename
        out_id_1 = "".join(
            random.choice(
                string.ascii_uppercase + string.digits + string.ascii_lowercase
            )
            for i in range(16)
        )

        outdir = os.path.abspath(f"{alignment_dir}/{out_id_1}")

        os.system(f"cp {contig1} {outdir}/1.tempcontig")
        os.system(f"cp {contig2} {outdir}/2.tempcontig")

        # Aligning the genomes
        os.system(
            f"cd {outdir}; mugsy "
            f"    --directory {outdir} "
            f"    --prefix out "
            f"    {outdir}/1.tempcontig "
            f"    {outdir}/2.tempcontig"
        )

        os.system(f"mv {outdir}/out.maf {final_name}")

        os.system(f"rm -rf {outdir}")

    return final_name


def calculate_length_bias(input_alignment, genome_1_file, genome_2_file, output_file):
    strain1, strain2 = input_alignment.rsplit("/", 1)[1].rsplit(".", 1)[0].split("_@_")
    g1size = sum(len(s) for s in SeqIO.parse(genome_1_file, "fasta"))
    g2size = sum(len(s) for s in SeqIO.parse(genome_2_file, "fasta"))
    (
        init_div,
        alignment_size,
        observed_sum_sq_diff,
        low_percentile,
        high_percentile,
    ) = get_transfer_measurement(input_alignment)

    edge = "\t".join(
        [
            strain1,
            strain2,
            str(init_div),
            str(alignment_size),
            str(g1size),
            str(g2size),
            str(observed_sum_sq_diff),
            str(low_percentile),
            str(high_percentile),
        ]
    )

    with open(output_file, "w") as outfile:
        outfile.write(edge + "\n")


def get_transfer_measurement(alignment: str, min_block_size=0, filtering_window=1000):
    # Initializes local variables
    all_blocks, prefilter_total_len = get_concatenated_alignment(alignment)

    # Filter alignment to split into subblocks at any point where there are at least 2 gaps
    filtered_blocks = [
        i
        for prefilter_s1, prefilter_s2 in all_blocks
        for i in filter_block(prefilter_s1, prefilter_s2)
    ]

    def blocks_alignment_div(filtered_blocks: Iterable[Tuple[str, str]]):
        s1temp, s2temp = zip(
            *(block for block in filtered_blocks if len(block[0]) > min_block_size)
        )

        # Assumes that each alignment block adds another divergence
        s1 = "1".join(s1temp)
        s2 = "0".join(s2temp)
        alignment_size = len(s1)
        init_div_count = naive_div_count(s1, s2)
        init_div = init_div_count * 1.0 / alignment_size

        return s1, s2, alignment_size, init_div

    _, _, _, init_div_raw = blocks_alignment_div(filtered_blocks)

    # Second filtering step by divergence
    final_filtered = [
        i
        for s1, s2 in filtered_blocks
        for i in filter_block_by_divergence(
            s1, s2, init_div_raw, winlen=filtering_window
        )
    ]
    s1, s2, alignment_size, init_div = blocks_alignment_div(final_filtered)

    if init_div > 0:
        initial = id_var_window_counts(s1, s2)
        initial_cumulative = get_cumulative_window_spectrum(initial, alignment_size)
        null_expect = single_param_null_model(
            np.arange(0, len(initial_cumulative)), init_div
        )
        observed_sum_sq_diff = np.sum(
            np.square(np.subtract(initial_cumulative, null_expect))
        )

        # Given a distribution of identical windows, bootstraps to find
        # length bias (SSD) confidence interval
        ssds = []
        for t in range(0, 200):
            initial_boot = np.random.choice(initial, len(initial), replace=True)
            initial_cumulative_boot = get_cumulative_window_spectrum(
                initial_boot, alignment_size  # type: ignore
            )
            ssd_boot = np.sum(
                np.square(np.subtract(initial_cumulative_boot, null_expect))
            )
            ssds.append(ssd_boot)
        low_percentile = np.percentile(ssds, 0.5)
        high_percentile = np.percentile(ssds, 99.5)
    else:
        observed_sum_sq_diff = np.nan
        low_percentile = np.nan
        high_percentile = np.nan

    return (
        init_div,
        alignment_size,
        observed_sum_sq_diff,
        low_percentile,
        high_percentile,
    )


def parse_alignment_file(alignment, min_block_size=1000, filtering_window=1000):
    # Initializes local variables
    all_blocks, prefilter_total_len = get_concatenated_alignment(alignment)

    # Filter alignment to split into subblocks at any point where there are at least 2 gaps
    filtered_blocks = [
        block
        for block in (
            filter_block(prefilter_s1, prefilter_s2)
            for prefilter_s1, prefilter_s2 in all_blocks
        )
        if len(block[0]) > min_block_size
    ]
    s1temp, s2temp = zip(*filtered_blocks)

    # Assumes that each alignment block adds another divergence
    Concat_S1 = "1".join(s1temp)
    Concat_S2 = "0".join(s2temp)
    alignment_size = len(Concat_S1)
    init_div_count = naive_div_count(Concat_S1, Concat_S2)
    init_div = init_div_count * 1.0 / alignment_size

    # Second filtering step by divergence
    final_filtered = []
    for s1, s2 in filtered_blocks:
        final_filtered += filter_block_by_divergence(
            s1, s2, init_div, winlen=filtering_window
        )
    filtered_blocks = [
        block for block in final_filtered if len(block[0]) > min_block_size
    ]
    s1temp, s2temp = zip(*filtered_blocks)

    # Assumes that each alignment block adds another divergence
    Concat_S1 = "1".join(s1temp)
    Concat_S2 = "0".join(s2temp)
    alignment_size = len(Concat_S1)
    init_div_count = naive_div_count(Concat_S1, Concat_S2)
    init_div = (init_div_count * 1.0) / alignment_size

    initial = id_var_window_counts(Concat_S1, Concat_S2)
    initial_cumulative = get_cumulative_window_spectrum(initial, alignment_size)

    return (initial_cumulative, init_div)


def get_cumulative_window_spectrum(idw: List[int], genome_size: int):
    """
    Gets the X and Y coordinates of the identical window spectrum
    i.e., the fraction of the genome belonging to identical windows
    above a certain size
    """
    obs_frac_counts = np.zeros(genome_size)
    norm = np.sum(idw)
    windows = Counter(idw)
    for wsize, count in windows.items():
        obs_frac_counts[wsize] = count * wsize * 1.0 / norm
    return 1.0 - np.cumsum(obs_frac_counts)


def get_concatenated_alignment(alignment: str):
    """
    This creates a list of tuples that constitute a concatenated alignment.
    Every entry in the list is a tuple that corresponds to an alignment block.
    """
    seqs: List[Tuple[str, str]] = []
    total_len = 0
    with open(alignment) as fi:
        """
        Parser assumes a maf format where every alignment block begins with a
        statement of how many sequences are in that block, indicated by
        "mult=." Also assumes that the order of sequences in each block is
        the same.
        """
        for line in fi:
            if "mult=2" in line:
                seq_line_1 = next(fi)
                block_1 = seq_line_1.split()[-1].strip()
                total_len += len(block_1)
                seq_line_2 = next(fi)
                block_2 = seq_line_2.split()[-1].replace("\n", "")
                seqs.append((block_1, block_2))
    return seqs, total_len


def id_var_window_counts(sequence_1, sequence_2):
    """
    This method takes two aligned sequences (strings) and returns the
    lengths of all identical windows between those sequences.
    """
    if sequence_1 == sequence_2:
        return [len(sequence_1)]
    a1 = np.array(list(sequence_1))
    a2 = np.array(list(sequence_2))
    mutated_positions = np.where(a1 != a2)[0]
    return (
        np.ediff1d(
            mutated_positions,
            to_begin=mutated_positions[0] + 1,
            to_end=len(sequence_1) - mutated_positions[-1],
        )
        - 1
    )


def naive_div_count(sequence_1, sequence_2):
    """
    Given two aligned strings, returns the number of differences between them.
    """
    if sequence_1 == sequence_2:
        return 0
    a1 = np.array(list(sequence_1))
    a2 = np.array(list(sequence_2))
    return len(np.where(a1 != a2)[0])


def filter_block(sequence_1, sequence_2):
    return [
        block
        for block in get_filtered_subblocks(
            sequence_1,
            sequence_2,
            filter_string(sequence_1) + filter_string(sequence_2),
        )
        if block[0] != ""
    ]


def filter_string(S: str):
    """
    >>> filter_string("ATG-CAT--GC--ATGC")
    [(9, 7), (13, 11)]
    """
    begin = 0
    filter_intervals: List[Tuple[int, int]] = []
    for base, count in ((label, sum(1 for _ in group)) for label, group in groupby(S)):
        if base == "-" and count >= 2:
            filter_intervals.append((begin + count, begin))
        begin += count
    return filter_intervals


def filter_block_by_divergence(sequence_1, sequence_2, init_div, winlen=1000):
    """
    Filters two sequences from an alignment block to remove regions
    that are significantly more diverged than expected
    """
    if sequence_1 == sequence_2:
        return [(sequence_1, sequence_2)]

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

    return [
        block
        for block in get_filtered_subblocks(sequence_1, sequence_2, removal_positions)
        if block[0] != ""
    ]


def get_filtered_subblocks(sequence_1, sequence_2, positions_to_remove):
    """
    Helper method that splits a string into substrings when given a list of
    start and end positions to remove
    """
    if len(positions_to_remove) == 0:
        return [(sequence_1, sequence_2)]
    final_blocks = []
    initial_start = 0
    for end, start in sorted(
        merge_intervals(sorted(positions_to_remove, reverse=True)), key=lambda x: x[1]
    ):
        yield sequence_1[initial_start:start], sequence_2[initial_start:start]
        initial_start = end
    yield sequence_1[initial_start:], sequence_2[initial_start:]


def single_param_null_model(w, div):
    """
    The simple single parameter null model that describes
    the window spectrum under an assumption of only mutation
    and no transfer
    """
    return np.exp(-div * w) * (div * w + 1)


def merge_intervals(intervals: List[Tuple[int, int]]):
    """
    >>> merge_intervals(sorted([(39, 37), (67, 65)], reverse=True))
    [(67, 65), (39, 37)]
    >>> merge_intervals(sorted([(3, 1), (4, 1)], reverse=True))
    [(4, 1)]
    """
    all_intervals = []
    current_interval = intervals[0]
    for j, interval in enumerate(intervals[1:]):
        end, start = interval
        if current_interval[1] <= end:
            current_interval = (current_interval[0], start)
        else:
            all_intervals.append(current_interval)
            current_interval = interval
    if current_interval not in all_intervals:
        all_intervals.append(current_interval)
    return all_intervals


def concat_length_bias_files(outfile: str, length_bias_files: List[str]):
    header = [
        "Strain 1",
        "Strain 2",
        "Initial divergence",
        "Alignment size",
        "Genome 1 size",
        "Genome 2 size",
        "Observed SSD",
        "SSD 95 CI low",
        "SSD 95 CI high",
    ]
    with open(outfile, "w") as fo:
        fo.write("\t".join(header) + "\n")
        for f in length_bias_files:
            with open(f) as fi:
                fo.write(fi.read())


if __name__ == "__main__":
    g1 = "../../test/M1612_contigs.fasta"
    g2 = "../../test/M1613_contigs.fasta"

    import random

    seed = random.randint(1, int(1e9))
    g1m = rename_for_mugsy(g1)
    g2m = rename_for_mugsy(g2)

    alignment_dir = "proc"
    maf = align_genomes(g1m, g2m, alignment_dir, seed)

    output = alignment_dir + "{g1}_@_{g2}.length_bias.txt"
    calculate_length_bias(maf, g1m, g2m, output)
