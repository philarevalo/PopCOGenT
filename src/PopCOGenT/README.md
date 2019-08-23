# PopCOGenT
Measures recent horizontal gene transfer between pairs of genomes using the length bias of the distribution of identical genome regions and produces population predictions.

# Important warnings

* PopCOGenT produces a lot of pairwise genome alignments. These are on the order of 10MB in size, so please use the `discard_alignments` option if you wish to save space.
* The config file has options for running on slurm. These are not functional (23 August 2019) and should be left as empty strings. If you want to run PopCOGenT on a HPCC, we recommend that you run PopCOGenT.sh as a batch script in your job management system of choice. Make sure that you request a node that has resources that match what you put in the config file (i.e., if you set the `-n` option to 4, request a node with 4 cores).

# Input
A directory of genomes to cluster in fasta format.

# Output

`*_length_bias.txt` contains the raw calculation of length bias between all genomes.
`*.graphml` contains the unclustered length bias network in graphml format for visualization.
`*.cluster.tab.txt` contains final population assignments.

# Setup
Set all relevant variables in `config.sh`. Variable descriptions are provided as comments in `config.sh` and are also reproduced below.

* `base_name`: Base name for output files. Just a prefix to identify your outputs.
* `final_output_dir`: Output directory for the final output files.
* `mugsy_path`: Path to mugsy. Please provide absolute path.
* `mugsy_env`: Path to `mugsyenv.sh`. Please provide absolute path.
* `infomap_path`: Path to Infomap executable. Please provide absolute path.
* `genome_dir`: Path to genome files.
* `genome_ext`: Genome file filename extension.
* `alignment_dir`: Directory for output alignments. Must provide absolute path.
* `num_threads`: Are you running on a single machine? Please specify the number of threads to run. This can, at maximum, be the number of logical cores your machine has.
* `single_cell`: Are your genomes single-cell genomes? If so, this should equal --single_cell. Otherwise leave as ''
* `discard_alignments`: Whether to discard alignments after length bias is calculated. Alignment files can be 10MB each and thus a run on 100 genomes can take up on the order of 50 GB of space if alignment files are not discarded. Please set as either `True` or `False`.

# Execution

Once all variables are set in `config.sh`, you can run PopCOGenT using the run script `PopCOGenT.sh`.
