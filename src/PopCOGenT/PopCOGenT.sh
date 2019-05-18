#!/bin/bash

configfile=./config.sh
source ${configfile}
source activate PopCOGenT
source ${mugsy_env}

if [ "${slurm_str}" = "" ]
	then
		python get_alignment_and_length_bias.py --genome_dir ${genome_dir} --genome_ext ${genome_ext} --alignment_dir ${alignment_dir} --mugsy_path ${mugsy_path} --mugsy_env ${mugsy_env} --base_name ${base_name} --final_output_dir ${final_output_dir} --num_threads ${num_threads} ${keep_alignments}
		python cluster.py --base_name ${base_name} --length_bias_file ${final_output_dir}/${base_name}.length_bias.txt --output_directory ${final_output_dir} --infomap_path ${infomap_path} ${single_cell}
fi
