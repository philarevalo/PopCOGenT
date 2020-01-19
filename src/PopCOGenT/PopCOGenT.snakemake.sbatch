#!/bin/bash

#SBATCH -n 1

source activate PopCOGenT
source /home/parevalo/apps/mugsy_trunk/mugsyenv.sh
snakemake --cluster-config cluster.yml --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem={cluster.mem}" --jobname {rulename}.{jobid} --jobs 250

