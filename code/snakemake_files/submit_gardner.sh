#!/bin/bash
snakemake -j $1  --cluster-config gardner.yaml --js gardner_script.sh --cluster "qsub -l walltime={cluster.time} -l nodes={cluster.nodes}:ppn={cluster.mincpus} -l mem={cluster.mem}" --rerun-incomplete -k
