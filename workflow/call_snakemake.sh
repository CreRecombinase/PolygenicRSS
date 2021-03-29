#!/bin/bash
. /gpfs/data/xhe-lab/nwk/miniconda3/etc/profile.d/conda.sh
module purge
conda activate snakemake
export SINGULARITY_CACHEDIR=/gpfs/data/xhe-lab/software/containers/
export SINGULARITY_DOCKER_USERNAME=crerecombinase
export SINGULARITY_DOCKER_PASSWORD=4151bbc097dd247eb8df2951c7d3deb19a7a2ea5
snakemake --profile ../config/gardner "$@" #|| ls -rt error/* | tail -1 | xargs cat
