# Path to genome files
genome_directory= # Insert path to genome files

# Extension of genome files
genome_extension= # insert file extension for genome files in fasta format

# Base name for output files
base_name= # Insert base filename identifier here

# Single cell genome flag
single_cell=""

# Output directory for final length bias file
final_output_dir=./output/results

# Directories for alignments and runscripts
temp_script_dir=./output/${base_name}_scripts/
slurm_output_dir=./output/${base_name}_ssd_slurm/
alignment_dir=./output/${base_name}_ssd_align/

# SLURM parameters for jobs
cluster_username= # insert cluster username here
partition_name= #insert name of partitions here
max_jobs=100
max_submit=50

# Path to mugsy
mugsy_path=/home/parevalo/apps/mugsy_trunk/mugsy
path_to_utilities=./
