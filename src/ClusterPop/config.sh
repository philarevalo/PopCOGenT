# Path to genome files
genome_directory=/nobackup1b/users/parevalo/genomes_for_clustering/inputs/SalTyphimurium/

# Extension of genome files
genome_extension=.fna

# Base name for output files
base_name=SalTyphimurium

# Single cell genome flag
single_cell=""

# Output directory for final length bias file
final_output_dir=/home/parevalo/scratch/HGT_Clustering/results/

# Directories for alignments and runscripts
temp_script_dir=/home/parevalo/scratch/genomes_for_clustering/${base_name}_scripts/
slurm_output_dir=/home/parevalo/scratch/genomes_for_clustering/${base_name}_ssd_slurm/
alignment_dir=/home/parevalo/scratch/genomes_for_clustering/${base_name}_ssd_align/

# SLURM parameters for jobs
cluster_username=parevalo
partition_name=sched_mit_chisholm,sched_mit_hill,newnodes
max_jobs=100
max_submit=50

# Path to mugsy
mugsy_path=/home/parevalo/apps/mugsy_trunk/mugsy
path_to_utilities=../utils/
