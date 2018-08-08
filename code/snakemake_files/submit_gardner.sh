#!/bin/bash
snakemake --verbose -j $1  --cluster-config gardner.yaml --cluster "qsub -V -l walltime={cluster.time} -l nodes={cluster.nodes}:ppn={cluster.cpuspertask} -l mem={cluster.mem}" --rerun-incomplete
