# ClusterPop
Measures recent horizontal gene transfer between pairs of genomes using the ClusterPop methodology and produces population predictions.

# Requirements
As written, this code requires use of a computational cluster environment with SLURM. Infomap is also required and not provided in this repository. A list of required python 3 packages is available in ClusterPop.yml. These should be installed using miniconda in an environment called ClusterPop.

# Input
A directory with genomes in fasta format.

# Output

`*_cat_ssd.txt` contains the raw calculation of length bias between all genomes.
`*.graphml` contains the unclustered length bias network in graphml format for visualization.
`*.cluster.tab.txt` contains final population assignments.

The base filename for all output files is set in `config.sh`.

# Usage
Set all relevant variables in `config.sh` and then run.

`sbatch run_ClusterPop.sbatch`