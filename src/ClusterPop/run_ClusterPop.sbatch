#!/bin/bash
source activate HGT_cluster
source ${mugsy_path}/mugsyenv.sh

if [ ! -d "$temp_script_dir" ]; then
     mkdir ${temp_script_dir}
fi

if [ ! -d "$slurm_output_dir" ]; then
     mkdir ${slurm_output_dir}
fi

if [ ! -d "$alignment_dir" ]; then
     mkdir ${alignment_dir}
fi

if [ ! -d "$final_output_dir" ]; then
     mkdir ${final_output_dir}
fi

time python ${path_to_utilities}/align_and_calc_length_bias.py --base_name ${base_name} --genome_dir ${genome_directory} --genome_ext ${genome_extension} --temp_script_dir ${temp_script_dir} --slurm_output_dir ${slurm_output_dir} --alignment_dir ${alignment_dir} --cluster_username ${cluster_username} --partition_name ${partition_name} --mugsy_path ${mugsy_path} --utility_path ${path_to_utilities} --max_jobs ${max_jobs} --max_submit ${max_submit} --final_output_dir ${final_output_dir}
time python ${path_to_utilities}/cluster.py -s ${final_output_dir}/${base_name}_cat_ssd.txt -o ${base_name} -d ${final_output_dir} ${single_cell}