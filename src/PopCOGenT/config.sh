# Base name for final output file
base_name='pcc'

# Path to mugsy and mugsyenv.sh. Please provide absolute path
mugsy_path=/home/parevalo/apps/mugsy_trunk/mugsy
mugsy_env=/home/parevalo/apps/mugsy_trunk/mugsyenv.sh

# Path to infomap. Please provide absolute path.
infomap_path=/nobackup1/parevalo/Infomap/Infomap

# Path to genome files
genome_dir=/nobackup1/parevalo/genomes_for_clustering/inputs/Pcc/clean/ #../../test/ #Insert path to genome files

# Are the genomes single-cell? If so, this should equal --single_cell
single_cell='--single_cell'

# Extension of genome files
genome_ext=.fasta # insert file extension for genome files in fasta format

# Directory for output alignments. Must provide absolute path.
alignment_dir=/home/parevalo/testing/${base_name}_ssd_align/

# Output directory for final length bias file
final_output_dir=./


# Are you running on a single machine? Please specify the number of threads
num_threads=4

# Are you using a slurm environment? Then this should equal --slurm, otherwise, leave as empty quotes.
slurm_str=''

# If using slurm, please specify the output directory for the runscripts and source scripts. Absolute paths required.
script_dir=''
source_path=''



source activate HGT_cluster
source ${mugsy_env}

if [ "${slurm_str}" = "" ]
	then
		python get_alignment_and_length_bias.py --genome_dir ${genome_dir} --genome_ext ${genome_ext} --alignment_dir ${alignment_dir} --mugsy_path ${mugsy_path} --mugsy_env ${mugsy_env} --base_name ${base_name} --final_output_dir ${final_output_dir} --num_threads ${num_threads}
		python cluster.py --base_name ${base_name} --length_bias_file ${final_output_dir}/${base_name}.length_bias.txt --output_directory ./infomap_out --infomap_path ${infomap_path} ${single_cell}
fi


