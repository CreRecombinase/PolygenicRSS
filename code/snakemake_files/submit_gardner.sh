#!/bin/bash

#--output={cluster.output} --error={cluster.error}
#--job-name={cluster.jobname}
snakemake --verbose -j $1  --cluster-config gardner.yaml --cluster "qsub -V -l walltime={cluster.time} -l nodes={cluster.nodes}:ppn={cluster.cpus} -l mem={cluster.mempercpu}" 
