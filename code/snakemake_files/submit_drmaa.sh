#!/bin/bash
snakemake --verbose -j $1  --cluster-config gardner.yaml --drmaa " -w -V -l walltime={cluster.time} -l nodes={cluster.nodes}:ppn={cluster.mincpus} -l mem={cluster.mem}" --rerun-incomplete -k --drmaa-log-dir=./drmaa_results/
