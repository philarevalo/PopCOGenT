#!/bin/bash

source activate HGT_cluster

module load engaging/mafft/7.245-with-extensions 

snakemake --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem={cluster.mem} --time={cluster.time}" --jobname {rulename}.{jobid} --jobs 150 --keep-going --local-cores 3
